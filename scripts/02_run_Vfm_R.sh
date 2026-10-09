#!/usr/bin/env bash

#SBATCH --job-name=Roux_Vfm_plot
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/Roux_Vfm_plot_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/Roux_Vfm_plot_%j.err

set -euo pipefail


# ============================================================
# VFM MANUSCRIPT FIGURES AND TABLES
#
# Runs:
#   scripts/02_Vfm_manuscript_outputs.R
#
# Creates the R environment automatically if it is absent.
# ============================================================


# ============================================================
# PATHS
# ============================================================

ROOT="/scratch/al98750/Roux"

PROJECT="${HOME}/Rb_Solanimycin"

ENV="${ROOT}/envs/vfm_R"

RSCRIPT="${PROJECT}/scripts/02_Vfm_manuscript_outputs.R"

OUTPUT="${ROOT}/06_QS_Vfm/07_manuscript"


# ============================================================
# JOB INFORMATION
# ============================================================

echo
echo "============================================================"
echo " VFM MANUSCRIPT FIGURE/TABLE ANALYSIS"
echo "============================================================"
echo

echo "Job ID:"
echo "  ${SLURM_JOB_ID:-interactive}"

echo
echo "Node:"
echo "  $(hostname)"

echo
echo "Started:"
echo "  $(date)"

echo


# ============================================================
# INITIALIZE CONDA
# ============================================================

# Avoid accidentally using the system/module R installation.
module unload R >/dev/null 2>&1 || true

source "$(conda info --base)/etc/profile.d/conda.sh"


# ============================================================
# CREATE R ENVIRONMENT IF ABSENT
# ============================================================

if [[ ! -x "${ENV}/bin/Rscript" ]]; then

    echo "============================================================"
    echo " R ENVIRONMENT NOT FOUND"
    echo " Creating:"
    echo "   ${ENV}"
    echo "============================================================"
    echo

    conda create \
        -p "${ENV}" \
        -c conda-forge \
        r-base=4.4 \
        r-ggplot2 \
        r-dplyr \
        r-tidyr \
        r-readr \
        r-forcats \
        r-patchwork \
        -y

else

    echo "R environment already exists:"
    echo "  ${ENV}"
    echo

fi


# ============================================================
# ACTIVATE R ENVIRONMENT
# ============================================================

conda activate "${ENV}"


echo "Using R:"
echo "  $(which R)"

echo
echo "Using Rscript:"
echo "  $(which Rscript)"

echo
Rscript --version
echo


# ============================================================
# CHECK REQUIRED R SCRIPT
# ============================================================

if [[ ! -s "${RSCRIPT}" ]]; then

    echo "ERROR:"
    echo "R script was not found:"
    echo "  ${RSCRIPT}"

    exit 1

fi


# ============================================================
# CHECK REQUIRED R PACKAGES
# ============================================================

echo "============================================================"
echo " Checking required R packages"
echo "============================================================"
echo

Rscript -e '
packages <- c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "readr",
    "forcats",
    "patchwork"
)

missing <- packages[
    !sapply(
        packages,
        requireNamespace,
        quietly = TRUE
    )
]

if (length(missing) > 0) {

    stop(
        paste(
            "Missing R packages:",
            paste(missing, collapse = ", ")
        )
    )

}

cat("All required R packages are available.\n")
'


# ============================================================
# RUN MANUSCRIPT ANALYSIS
# ============================================================

echo
echo "============================================================"
echo " Running Vfm manuscript analysis"
echo "============================================================"
echo

echo "R script:"
echo "  ${RSCRIPT}"

echo
echo "Output directory:"
echo "  ${OUTPUT}"

echo


cd "${PROJECT}"


Rscript "${RSCRIPT}"


# ============================================================
# SHOW CREATED OUTPUTS
# ============================================================

echo
echo "============================================================"
echo " GENERATED FILES"
echo "============================================================"
echo

if [[ -d "${OUTPUT}" ]]; then

    find "${OUTPUT}" \
        -maxdepth 2 \
        -type f \
        -printf "%p\n" \
        | sort

else

    echo "WARNING:"
    echo "Output directory was not created:"
    echo "  ${OUTPUT}"

fi


# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================================"
echo " VFM MANUSCRIPT ANALYSIS COMPLETE"
echo "============================================================"
echo

echo "Finished:"
echo "  $(date)"

echo