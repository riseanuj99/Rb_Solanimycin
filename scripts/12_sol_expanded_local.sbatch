#!/bin/bash
#SBATCH --job-name=sol_local_discovery
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=sol_local_discovery_%j.out
#SBATCH --error=sol_local_discovery_%j.err
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /scratch/al98750/Roux/envs/sol_hgt
cd "$HOME/Rb_Solanimycin/scripts"
bash 12_sol_expanded_local.sh analyze
