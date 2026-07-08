import yaml
from pathlib import Path
import pandas as pd
import re
import numpy as np
from typing import Optional, Sequence, Callable
import copy
from tqdm import tqdm 

import optuna as opt
from optuna.trial import Trial, FrozenTrial
from optuna.study import Study

from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import VarianceThreshold, SelectKBest
from sklearn.utils.class_weight import compute_sample_weight
import sklearn.metrics as metrics
import xgboost as xgb
from xgboost import DMatrix, QuantileDMatrix
import shap

DROP = {
    "Niche": ['Ear', 'Sputum', 'Plant', 'Skin', 'Throat', 'Upper.Respiratory'],
    "Cluster": []
    }
FEATURES = Path("")
LABELS = Path("")
SPLITS = Path("")
HYPOTHESIS_FILE = Path("")
OUTPUT = Path("")

# --------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------
def load_data(FEATURES, LABELS, SPLITS):
    """
    Load files once - features are huge so reading is slow
    """
    X = pd.read_parquet(FEATURES).sort_index()
    labels = pd.read_parquet(LABELS).sort_index()
    splits = pd.read_parquet(SPLITS)
    return {
        "Features": X,
        "Labels": labels,
        "Splits": splits
    }

def prep_data(X, labels, splits, drops):
    """
    Make the train and val splits based on the `drop` instructions.
    """
    features = X.columns.to_numpy()

    if drops is not None:
        labels = labels.loc[X.index]
        for l, d in drops.items():
            labels = labels[~labels[l].isin(d)]
        X = X.loc[labels.index] 
        splits = splits.loc[labels.index] 

    # Make splits  
    train_mask = splits["split"] == "Train"
    val_mask = splits["split"] == "Val"
    test_mask = splits["split"] == "Test"

    X_train = X[train_mask].to_numpy()
    y_train = labels[train_mask]
    X_val = X[val_mask].to_numpy()
    y_val = labels[val_mask]
    X_test = X[test_mask].to_numpy()
    y_test = labels[test_mask]
    return {
        "Raw": (X, labels),
        "Train": (X_train, y_train), 
        "Val": (X_val, y_val),
        "Test": (X_test, y_test),
        "Features": features
        }

def read_yaml(path: Path):
    path = path.with_suffix(".yaml")
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)

HYPOTHESIS = read_yaml(HYPOTHESIS_FILE)

NUM_CLASSES = len(set(HYPOTHESIS.values()))

# Training device
GPU_ID = -1

# Model parameters 
MPARAM = {
    "booster": "gbtree",
    "device": "cpu" if GPU_ID < 0 else f"cuda:{GPU_ID}",
    "eval_metric": "mlogloss",
    "objective": "multi:softmax",
    "num_class": NUM_CLASSES,
    "tree_method": "hist",
    "verbosity": 0,
    # ---
    "base_score": None,
    "colsample_bylevel": None,
    "colsample_bynode": None,
    "grow_policy": None,
    "interaction_constraints": None,
    "max_bin": None,
    "max_cat_threshold": None,
    "max_cat_to_onehot": None,
    "max_delta_step": None,
    "max_leaves": None,
    "monotone_constraints": None,
    "multi_strategy": None,
    "num_parallel_tree": None,
    "random_state": None,
    "reg_alpha": None,
    "sampling_method": None,
    "validate_parameters": None,
}

# Dmatrix parameters 
DPARAM = {
    "max_bin": None,
    "n_jobs": -1,
}

# Training parameters 
TPARAM = {
    "num_boost_round": 100,
    "evals": [],
    "early_stopping_rounds": None,
    "evals_result": {},
    "obj": None,
    "custom_metric": None,
    "verbose_eval": True,
    "callbacks": None,
}

# Default hyperparameters
HPARAM = {
    "optimize":0,
    "learning_rate":1e-2,
    "max_depth":3,
    "min_child_weight":1,
    "gamma":1e-8,
    "subsample":0.5,
    "colsample_bytree":0.5,
    "reg_lambda":1e-8,
    "alpha":1e-8,
    "scale_pos_weight":1e-3,
}

# --------------------------------------------------------------------
# Data processing 
# --------------------------------------------------------------------
def label_mapping(HYPOTHESIS, labels):
    remap = {l: HYPOTHESIS.get(l, l) for l in labels["Niche"]}
    classes = list(dict.fromkeys(remap.values()))
    le = LabelEncoder().fit(classes)
    encode = lambda c: int(le.transform([c])[0])
    encoder = {l: encode(c) for l, c in remap.items()}
    decoder = {encode(c): c for c in classes}
    return encoder, decoder

def encode_labels(encoder, y):
    return y["Niche"].map(encoder).to_numpy()

def decode_labels(decoder, codes):
    return pd.Series(codes).map(decoder).to_numpy()

def best_features(X, y, encoder, k: int = 20000):
    """Reduce features to top k ANOVA F-scores."""
    # Compute non-constant features
    non_constant = VarianceThreshold().fit(X, y).get_support()
    X = X[:, non_constant]

    # Compute high ANOVA F-value features
    k = min(k, X.shape[1])
    high_anova_f = SelectKBest(k=k).fit(X, y)
    high_anova_f = high_anova_f.get_support()

    # Intersect features
    nc_subset = non_constant[non_constant]
    nc2anova_map = nc_subset * high_anova_f
    mask = non_constant
    mask[np.where(non_constant)[0]] = nc2anova_map
    return mask

def create_dmatrix(
    X: np.ndarray,
    y: Optional[np.ndarray],
    feature_names: Optional[Sequence[str]] = None,
    weight: Optional[np.ndarray] = None,
    max_bin: Optional[int] = None,
    n_jobs: int = -1,
) -> DMatrix:
    param = {
        "data": X,
        "label": y,
        "feature_names": feature_names,
        "weight": weight,
        "nthread": n_jobs,
    }
    try:
        return QuantileDMatrix(**param, max_bin=max_bin, ref=None)
    except TypeError:
        return DMatrix(**param)


# --------------------------------------------------------------------
# Hyperparameter Tuning
# --------------------------------------------------------------------
class Fold(object):
    def __init__(self, X: np.ndarray, y: np.ndarray):
        self.X, self.y = X, y
        self.mask = np.zeros(len(self.y), dtype=bool)

    def extend(self, idx: np.ndarray):
        """Extend fold subsample with samples at idx."""
        self.mask[idx] = True

    @property
    def train(self) -> xgb.DMatrix:
        X = self.X[~self.mask]
        y = self.y[~self.mask]
        w = compute_sample_weight(class_weight="balanced", y=y)
        return create_dmatrix(X, y, weight=w)

    @property
    def validate(self) -> xgb.DMatrix:
        X = self.X[self.mask]
        y = self.y[self.mask]
        return create_dmatrix(X, y)


def stratified_kfold(
    encoder, 
    X: np.ndarray, 
    y: pd.DataFrame, 
    splits: int = 5
    ) -> list[Fold]:
    """Split data into k folds with uniform class proportions."""
    features = best_features(X, y, encoder)
    X = X[:, features]
    y = encoder(y)
    folds = [Fold(X, y) for _ in range(splits)]

    for code in range(len(label_group.classes)):
        # Find instances of category in data
        idx = np.argwhere(y == code)
        # Evenly split instances randomly
        np.random.shuffle(idx)
        idxs = np.array_split(idx, splits)
        # Extend current folds with new category
        for fold, idx in zip(folds, idxs):
            fold.extend(idx)

    return folds


def print_trial(_: Study, trial: FrozenTrial):
    """Helper utility to view trial statistics."""
    click.echo(f"Trial {trial.number}:")
    click.echo(f"  Value: {trial.value}")
    click.echo(f"  Params: ")
    for key, value in trial.params.items():
        click.echo(f"    {key}: {value}")

def sample_hparam(trial: Trial):
    return {
        "learning_rate":trial.suggest_float(
            "learning_rate", 0.01, 0.3, log=True
        ),
        "max_depth":trial.suggest_int("max_depth", 3, 20),
        "min_child_weight":trial.suggest_int("min_child_weight", 1, 10),
        "gamma":trial.suggest_float("gamma", 1e-8, 1.0, log=True),
        "subsample":trial.suggest_float("subsample", 0.5, 1.0),
        "colsample_bytree":trial.suggest_float("colsample_bytree", 0.5, 1.0),
        "reg_lambda":trial.suggest_float("reg_lambda", 1e-8, 1.0, log=True),
        "alpha":trial.suggest_float("alpha", 1e-8, 1.0, log=True),
        "scale_pos_weight":trial.suggest_float(
            "scale_pos_weight", 1e-3, 1.0, log=True
        ),
    }

def optimize(
    encoder,
    hparam: dict,
    X: np.ndarray,
    y: pd.DataFrame,
    gpu_id: int,
    verbose: bool,
):
    """Optimize hyperparameters using cross-validation."""
    # Set verbosity
    log_kind = opt.logging.INFO if verbose else opt.logging.WARNING
    log_callback = [print_trial] if verbose else []
    opt.logging.set_verbosity(log_kind)

    # Build optimizer
    pruner = opt.pruners.MedianPruner(n_warmup_steps=10)
    sampler = opt.samplers.TPESampler(seed=42)
    study = opt.create_study(
        study_name="optuna-optimizer",
        direction="maximize",
        pruner=pruner,
        sampler=sampler,
    )

    # Initialize trials with default hparams and folds
    study.enqueue_trial(dict(**hparam))
    folds = stratified_kfold(encoder, X, y, splits=5)

    # Optimize hyperparameters
    study.optimize(
        lambda trial: objective(trial, label_group, folds, gpu_id),
        n_trials=hparam["optimize"],
        n_jobs=-1,
        callbacks=log_callback,
        show_progress_bar=verbose,
    )

    best = study.best_params | {"optimize": hparam["optimize"]}
    return best

# --------------------------------------------------------------------
# Training
# --------------------------------------------------------------------
def train( 
    X: np.ndarray, 
    features: np.ndarray,
    y: pd.DataFrame, 
    k: int, 
    encoder: Callable,
    hparam: dict
    ):
    # Tune hyperparameters if specified

    # Filter features 
    y = encode_labels(encoder, y)
    top_features = best_features(X, y, encoder, k)
    X = X[:, top_features]
    f = features[top_features].tolist()
    w = compute_sample_weight(class_weight="balanced", y=y)

    dtrain = create_dmatrix(X, y, f, w, **DPARAM)
    params = MPARAM | DPARAM | hparam

    clf = xgb.train(
            params,
            dtrain,
            xgb_model=None,
            **TPARAM,
        )
    return top_features, clf


# --------------------------------------------------------------------
# Evaluation
# --------------------------------------------------------------------
def micro_f1_score(true, pred) -> float:
    return metrics.f1_score(true, pred, average="weighted")


def macro_f1_score(true, pred) -> float:
    return metrics.f1_score(true, pred, average="macro")


def accuracy_score(true, pred) -> float:
    return metrics.accuracy_score(true, pred)


def balanced_accuracy_score(true, pred) -> float:
    return metrics.balanced_accuracy_score(true, pred, adjusted=True)

def eval(
    clf,
    X: np.ndarray, 
    y: pd.DataFrame,
    top_features: np.ndarray,
    encoder: Callable,
    *metrics: Callable
    ):
    X = X[:, top_features]
    y = encode_labels(encoder, y)
    f = clf.feature_names
    deval = create_dmatrix(X, y, f, **DPARAM)
    true = deval.get_label()
    pred = clf.predict(deval)
    metrics = map(lambda dist: dist(true, pred), metrics)
    return pred, metrics


# --------------------------------------------------------------------
# Analysis 
# --------------------------------------------------------------------
def filter_genomovars(labels: pd.DataFrame, threshold: int = 50, genomovar_col: str = "Cluster"): 
    """
    Get genomovars from labels and return them based on a sample threshold. 
    Groups with a tiny number of samples are unreliable for evaluation. 
    """
    all_genomovar = labels[genomovar_col].unique()
    return [g for g in all_genomovar if (labels[genomovar_col] == g).sum()]

def calculate_shap(
    model,
    encoder,
    X,
    y,
    features,
    feature_names,
    dparam
    ):
    # Set model and inputs 
    # trained_error = "Model is untrained - cannot calculate SHAP values"
    # assert not (model.clf is None or model.features is None), trained_error
    X = X[:, features]
    X = X.reshape(X.shape[0], X.shape[-1])
    y = encode_labels(encoder, y)
    f = feature_names
    dmatrix = create_dmatrix(X, y, f, **dparam)
    # Shap explanation
    explainer = shap.TreeExplainer(model)
    explanation = explainer(dmatrix)
    explanation.data = explanation.data.toarray()
    explanation.feature_names = f
    return explanation

def shap_to_df(
    model,
    name,
    explanation, 
    labels: pd.DataFrame, 
    save_path: Path
    ):
    """
    Saves the SHAP values of an evaluated model. SHAP values are saved 
    as a set csv file. 1 csv for every label the model classifies to. 
    So, if there are 10 classes and 1000 samples, there will be 10 csv
    files with 1000 samples worth of SHAP values each. 
    """
    remap = {l: HYPOTHESIS.get(l, l) for l in labels["Niche"]}
    classes = list(dict.fromkeys(remap.values()))
    genomes = labels.index
    n_labels = explanation.values.shape[2]
    all_frames = []
    for n in range(n_labels):
        col_names = explanation.feature_names
        frame = pd.DataFrame(explanation.values[:,:,n], columns=col_names)
        frame.index = genomes
        frame.to_csv(save_path / 'shap_features_{}_niche_{}.csv'.format(name, sorted(classes)[n].replace("/", "-")))
        all_frames.append(frame)
    return all_frames

def holdout_loop(
    hparam: dict,
    holdout_genomovars: list[str],
    output: str | Path
    ):
    """
    Trains a model for each genomovar provided in ``holdout_genomovars``. Each 
    model is trained blind to the reference genomovar so it can be evaluated on it 
    after. This allows us to determine if the models are able to generalize to unseen
    genomovars. After evaluation, SHAP values are calculated to determine feature 
    importance for different environment niches. 

    Args:
        hparam (dict): Model hyperparameters 
        holdout_genomovars (list[str]): List of genomovars to filter to 
        output (str | Path): Folder for output artifacts
    """
    output = Path(output) # str to Path if needed
    loaded_data = load_data(FEATURES, LABELS, SPLITS)
    # Find all genomovars across all labels 
    all_genomovar = loaded_data["Labels"]['Cluster'].unique().tolist()
    single_cluster_filter = copy.deepcopy(DROP)     # A filter that removes the target genomovar  
    holdout_filter = copy.deepcopy(DROP)            # A filter that removes all non-target genomovar 

    metrics_df = pd.DataFrame(
        columns=[
            "Micro F1",
            "Macro F1",
            "Accuracy",
            "Balanced Accuracy",
            "True Labels",
            "Predicted Labels",
            ]
        )
    for genomovar in tqdm(holdout_genomovars):
        # Load data splits for all data not in `genomovar`
        single_cluster_filter["Cluster"] = [genomovar]
        holdout_filter["Cluster"] = [alt_g for alt_g in all_genomovar if alt_g != genomovar]
        main_data = prep_data(loaded_data["Features"], loaded_data["Labels"], loaded_data["Splits"], single_cluster_filter)
        X_train, y_train = main_data["Train"]
        X_val, y_val = main_data["Val"]
        encoder, decoder = label_mapping(HYPOTHESIS, y_train)
        features = data_sets["Features"]

        # Train the current model
        top_features, model = train(
            X_train, 
            features,
            y_train, 
            k = 20000, 
            encoder = encoder, 
            hparam = hparam
        )

        model.save_model(f"{genomovar}_holdout.ubj")

        pred, (f1, mf1, acc, bacc) = eval(
            model,
            X_val, 
            y_val,
            top_features,
            encoder,
            micro_f1_score,
            macro_f1_score,
            accuracy_score,
            balanced_accuracy_score
            )

        true = encode_labels(encoder, y_val)
        true_labels = y_val
        pred_labels = decode_labels(decoder, pred)

        metrics_df.loc[genomovar + "_baseline"] = [
            f1,
            mf1,
            acc,
            bacc,
            str(list(true_labels)),
            str(list(pred_labels)),
        ]

        # Load the target genomovar 
        holdout_data = prep_data(loaded_data["Features"], loaded_data["Labels"], loaded_data["Splits"], holdout_filter)
        # Stack all samples for validation 
        X_holdout = np.vstack((holdout_data["Train"][0], holdout_data["Val"][0], holdout_data["Test"][0]))
        y_holdout = pd.concat([holdout_data["Train"][1], holdout_data["Val"][1], holdout_data["Test"][1]])

        if len(y_holdout) == 0:
            print(f"{genomovar} has no samples in Clade A - skipping")
            continue

        pred, (f1, mf1, acc, bacc) = eval(
            model,
            X_holdout, 
            y_holdout,
            top_features,
            encoder,
            micro_f1_score,
            macro_f1_score,
            accuracy_score,
            balanced_accuracy_score
            )

        true = encode_labels(encoder, y_holdout)
        true_labels = y_holdout
        pred_labels = decode_labels(decoder, pred)

        metrics_df.loc[genomovar + "_holdout"] = [
            f1,
            mf1,
            acc,
            bacc,
            str(list(true_labels)),
            str(list(pred_labels)),
        ]

        # Perform SHAP evaluation on all validation samples - holdout + baseline 
        X_full_eval = np.vstack((X_val, X_holdout))
        y_full_eval = pd.concat([y_val, y_holdout])
        shap_values = calculate_shap(
            model,
            encoder,
            X_full_eval,
            y_full_eval,
            top_features,
            features[top_features].tolist(),
            DPARAM
        )

        # Save SHAP values 
        shap_to_df(
            model,
            genomovar,
            shap_values, 
            y_full_eval, 
            output
            )

    # Save metrics across all holdouts 
    metrics_df.to_parquet(output / "holdout_stats.parquet")


if __name__ == "__main__":
    loaded_data = load_data(FEATURES, LABELS, SPLITS)
    data_sets = prep_data(loaded_data["Features"], loaded_data["Labels"], loaded_data["Splits"], DROP)
    X_train, y_train = data_sets["Train"]
    features = data_sets["Features"]
    encoder, decoder = label_mapping(HYPOTHESIS, y_train)

    if HPARAM["optimize"] > 0:
        hparam = optimize(
            encoder,
            HPARAM,
            X_train,
            y_train,
            GPU_ID,
            verbose=True,
        )
    else:
        hparam = HPARAM

    # top_features, model = train(
    #     X_train, 
    #     features,
    #     y_train, 
    #     k = 20000, 
    #     encoder = encoder, 
    #     hparam = hparam
    # )
    # model.save_model("test_xgboost.ubj")

    # X_test, y_test = data_sets["Val"]
    # pred, (f1, mf1, acc, bacc) = eval(
    #     model,
    #     X_test, 
    #     y_test,
    #     top_features,
    #     encoder,
    #     micro_f1_score,
    #     macro_f1_score,
    #     accuracy_score,
    #     balanced_accuracy_score
    #     )

    # print(f"Micro F1: {f1}")
    # print(f"Macro F1: {mf1}")
    # print(f"Accuracy: {acc}")
    # print(f"Balanced Accuracy: {bacc}\n")

    threshold=50
    holdout_genomovars = filter_genomovars(loaded_data["Labels"], threshold)
    print(f"Found {len(holdout_genomovars)} genomovars with > {threshold} samples.")

    holdout_loop(
        hparam,
        holdout_genomovars=holdout_genomovars,
        output=OUTPUT,
        )