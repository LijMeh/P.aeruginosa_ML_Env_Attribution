#!/bin/bash
# ============================================================================
# Process Panaroo Output
# ============================================================================
# Processes the output from Panaroo and prepares it for model input.
# ============================================================================

eval "$(micromamba shell hook --shell bash)"

micromamba activate Protect

script_path="00_data/P.aeruginosa_ML_Env_Attribution/04_model_input_prep/00_model_input_prep"
input_path="00_data/panarooOutput"
output_path="00_data/model_input"

if [ ! -d "${output_path}" ]; then
    mkdir -p "${output_path}"
fi

echo "Processing Panaroo output..."

python ${script_path}/00_processPanarooOutput.py ${input_path} ${output_path}

echo "Removing duplicates..."

python ${script_path}/00_removeDuplicates.py ${output_path}

