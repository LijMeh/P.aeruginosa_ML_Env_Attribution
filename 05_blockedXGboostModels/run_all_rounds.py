"""
Train and export models for all 6 rounds.
Assumes parquet files exist:
  - features.parquet
  - labels.parquet  
  - splits.parquet
"""
from pathlib import Path
import sys
from trainBlocks import holdout_loop, filter_genomovars, load_data, read_yaml

# Paths to input data
FEATURES = Path("features.parquet")
LABELS = Path("labels.parquet")
SPLITS = Path("splits.parquet")

# Output directory
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Default hyperparameters
HPARAM = {
    "optimize": 0,
    "learning_rate": 1e-2,
    "max_depth": 3,
    "min_child_weight": 1,
    "gamma": 1e-8,
    "subsample": 0.5,
    "colsample_bytree": 0.5,
    "reg_lambda": 1e-8,
    "alpha": 1e-8,
    "scale_pos_weight": 1e-3,
}

def run_round(round_num: int, hparam: dict):
    """Train and export for a single round."""
    round_name = f"round{round_num}"
    
    # Load hypothesis for this round
    hypothesis_file = Path(f"{round_name}.yaml")
    if not hypothesis_file.exists():
        print(f"ERROR: {hypothesis_file} not found")
        return False
    
    hypothesis = read_yaml(hypothesis_file)
    
    # Create output directory for this round
    round_output = OUTPUT_DIR / round_name
    round_output.mkdir(exist_ok=True)
    
    # Load data
    try:
        loaded_data = load_data(FEATURES, LABELS, SPLITS)
    except Exception as e:
        print(f"ERROR loading data: {e}")
        return False
    
    # Get genomovars to holdout
    threshold = 50
    holdout_genomovars = filter_genomovars(loaded_data["Labels"], threshold)
    print(f"{round_name}: Found {len(holdout_genomovars)} genomovars with > {threshold} samples")
    
    if not holdout_genomovars:
        print(f"ERROR: No genomovars found for {round_name}")
        return False
    
    # Train and export
    try:
        holdout_loop(
            hparam=hparam,
            holdout_genomovars=holdout_genomovars,
            output=round_output,
        )
        print(f"✓ {round_name} completed successfully")
        return True
    except Exception as e:
        print(f"✗ {round_name} failed: {e}")
        return False

if __name__ == "__main__":
    print("Training and exporting models for 6 rounds...\n")
    
    results = {}
    for round_num in range(1, 7):
        results[f"round{round_num}"] = run_round(round_num, HPARAM)
    
    print("\n" + "="*50)
    print("Summary:")
    for round_name, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"  {round_name}: {status}")
    
    all_passed = all(results.values())
    sys.exit(0 if all_passed else 1)
