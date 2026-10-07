#!/bin/bash
#SBATCH --job-name=HiVir_cov_plot
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/al98750/Pana_RNAseq/logs/HiVir_cov_plot_%j.out
#SBATCH --error=/scratch/al98750/Pana_RNAseq/logs/HiVir_cov_plot_%j.err
#SBATCH --mail-user=al98750@uga.edu
#SBATCH --mail-type=END,FAIL


# =====================================================================
# HiVir RNA-seq coverage plotting
#
# This script:
#
#   1. Loads R
#   2. Creates/uses a personal R package library
#   3. Installs ggplot2 if it is missing
#   4. Checks all required files
#   5. Runs 06_plot_HiVir_coverage.R
#   6. Checks that PNG and PDF were created
#
# =====================================================================


set -e


# =====================================================================
# 1. PATHS
# =====================================================================

PROJECT="/scratch/al98750/Pana_RNAseq"

ANALYSIS_DIR="$PROJECT/HiVir_promoter_analysis"

R_SCRIPT="$ANALYSIS_DIR/06_plot_HiVir_coverage.R"

INPUT_FILE="$ANALYSIS_DIR/coverage/04_genotype_mean_pepM_coverage.csv"

PLOT_DIR="$ANALYSIS_DIR/plots"

LOG_DIR="$PROJECT/logs"


mkdir -p "$LOG_DIR"

mkdir -p "$PLOT_DIR"



# =====================================================================
# 2. START INFORMATION
# =====================================================================

echo
echo "============================================================"
echo "HiVir RNA-seq coverage plotting"
echo "============================================================"
echo

echo "Start time: $(date)"

echo "Compute node: $(hostname)"

echo "Analysis directory:"
echo "$ANALYSIS_DIR"

echo



# =====================================================================
# 3. LOAD R
# =====================================================================

module purge

module load R


echo "R version:"

R --version | head -1


echo

echo "Rscript location:"

which Rscript


echo



# =====================================================================
# 4. PERSONAL R PACKAGE LIBRARY
# =====================================================================
#
# Packages installed here will remain available for future jobs.
#
# =====================================================================

R_LIB="$HOME/R/x86_64-pc-linux-gnu-library/4.5"


mkdir -p "$R_LIB"


export R_LIBS_USER="$R_LIB"


echo "Personal R library:"

echo "$R_LIBS_USER"

echo



# =====================================================================
# 5. CHECK REQUIRED INPUT FILES
# =====================================================================

if [ ! -d "$ANALYSIS_DIR" ]
then

    echo "ERROR:"
    echo "Analysis directory does not exist:"
    echo "$ANALYSIS_DIR"

    exit 1

fi


if [ ! -f "$R_SCRIPT" ]
then

    echo "ERROR:"
    echo "R script not found:"
    echo "$R_SCRIPT"

    exit 1

fi


if [ ! -f "$INPUT_FILE" ]
then

    echo "ERROR:"
    echo "Coverage input file not found:"
    echo "$INPUT_FILE"

    exit 1

fi



# =====================================================================
# 6. INSTALL / CHECK REQUIRED R PACKAGES
# =====================================================================
#
# ggplot2 dependencies are installed automatically.
#
# We DO NOT use dependencies=TRUE because that also installs many
# optional suggested packages that are not needed here.
#
# =====================================================================

echo
echo "============================================================"
echo "Checking R packages"
echo "============================================================"
echo


Rscript - <<'RSCRIPT'

user_lib <- Sys.getenv("R_LIBS_USER")


dir.create(
  user_lib,
  recursive = TRUE,
  showWarnings = FALSE
)


.libPaths(
  c(
    user_lib,
    .libPaths()
  )
)


cat(
  "\nR library paths:\n"
)

print(
  .libPaths()
)


required_packages <- c(
  "ggplot2"
)


for (pkg in required_packages) {

  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {

    cat(
      "\nInstalling ",
      pkg,
      "...\n",
      sep = ""
    )


    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      lib = user_lib,
      Ncpus = 2
    )

  } else {

    cat(
      "\nAlready installed: ",
      pkg,
      "\n",
      sep = ""
    )
  }
}



# ---------------------------------------------------------
# Verify installation
# ---------------------------------------------------------

for (pkg in required_packages) {

  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {

    stop(
      paste(
        "Package installation failed:",
        pkg
      )
    )
  }


  cat(
    pkg,
    " version: ",
    as.character(
      packageVersion(
        pkg
      )
    ),
    "\n",
    sep = ""
  )
}


cat(
  "\nAll required R packages are available.\n"
)

RSCRIPT



# =====================================================================
# 7. MOVE TO ANALYSIS DIRECTORY
# =====================================================================

cd "$ANALYSIS_DIR"


echo

echo "Current working directory:"

pwd

echo



# =====================================================================
# 8. FIX OUTPUT PATH IF OLD VERSION OF R SCRIPT IS PRESENT
# =====================================================================
#
# Earlier version contained:
#
#     ../plots/
#
# We want:
#
#     plots/
#
# because this job runs from HiVir_promoter_analysis/
#
# =====================================================================

sed -i \
    's#"../plots/#"plots/#g' \
    "$R_SCRIPT"



# =====================================================================
# 9. SHOW THE INPUT DATA SUMMARY
# =====================================================================

echo
echo "============================================================"
echo "Coverage input check"
echo "============================================================"
echo


echo "First five lines:"

head -5 "$INPUT_FILE"


echo

echo "Number of positions per genotype:"


cut -d',' -f1 "$INPUT_FILE" \
    | tail -n +2 \
    | sort \
    | uniq -c


echo



# =====================================================================
# 10. RUN THE R COVERAGE PLOT
# =====================================================================

echo
echo "============================================================"
echo "Running HiVir coverage plotting script"
echo "============================================================"
echo


Rscript "$R_SCRIPT"



# =====================================================================
# 11. CHECK OUTPUTS
# =====================================================================

PNG_FILE="$PLOT_DIR/HiVir_pepM_RNAseq_coverage.png"

PDF_FILE="$PLOT_DIR/HiVir_pepM_RNAseq_coverage.pdf"


echo
echo "============================================================"
echo "Checking plot outputs"
echo "============================================================"
echo


if [ -f "$PNG_FILE" ]
then

    echo "PNG created successfully:"

    ls -lh "$PNG_FILE"

else

    echo "ERROR:"
    echo "PNG file was not created."

    exit 1

fi


echo


if [ -f "$PDF_FILE" ]
then

    echo "PDF created successfully:"

    ls -lh "$PDF_FILE"

else

    echo "ERROR:"
    echo "PDF file was not created."

    exit 1

fi



# =====================================================================
# 12. FINAL MESSAGE
# =====================================================================

echo
echo "============================================================"
echo "HiVir coverage analysis finished successfully"
echo "============================================================"
echo

echo "Plots are located in:"

echo "$PLOT_DIR"

echo

echo "PNG:"

echo "$PNG_FILE"

echo

echo "PDF:"

echo "$PDF_FILE"

echo

echo "End time: $(date)"

echo "============================================================"