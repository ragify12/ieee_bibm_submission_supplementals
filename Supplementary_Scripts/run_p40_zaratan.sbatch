#!/bin/bash
#SBATCH --job-name=p40_screen
#SBATCH --partition=gpu              # <-- your Zaratan GPU partition
#SBATCH --gres=gpu:a100:1            # one A100 (same as the paper)
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=p40_screen_%j.log

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- 1. environment -------------------------------------------------------
# Activate the SAME ColabFold env you used for the paper. Replace these two
# lines with whatever you already load on Zaratan (module load / conda activate).
module load cuda                       # <-- adjust to your Zaratan module
source activate colabfold              # <-- your ColabFold conda env

# --- 2. build binder:p40 FASTAs ------------------------------------------
python leads_vs_p40.py build

# --- 3. predict the 5 binder+p40 complexes (multimer, MSA on) -------------
colabfold_batch \
    --num-models 5 \
    --num-recycle 3 \
    p40_inputs/ p40_out/

# --- 4. ipSAE on each prediction + final table ---------------------------
python leads_vs_p40.py collect

echo "DONE. Result: p40_selectivity.csv"
