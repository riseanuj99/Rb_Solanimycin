#!/usr/bin/env bash
#SBATCH --job-name=solHGT_trees
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/solHGT_trees_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/solHGT_trees_%j.err
set -euo pipefail
ROOT=/scratch/al98750/Roux
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ROOT/envs/sol_hgt"
for x in orthofinder diamond mafft trimal iqtree2 python; do
    command -v "$x" >/dev/null || { echo "Missing: $x" >&2; exit 1; }
done
python "$HOME/Rb_Solanimycin/scripts/03_sol_HGT.py" phylogeny --root "$ROOT" --threads "${SLURM_CPUS_PER_TASK:-8}"
