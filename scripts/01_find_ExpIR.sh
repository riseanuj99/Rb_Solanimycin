#!/usr/bin/env bash

#SBATCH --job-name=Roux_ExpIR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/Roux_ExpIR_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/Roux_ExpIR_%j.err

set -euo pipefail
shopt -s nullglob


# ============================================================
# OBJECTIVE 1
# AHL QUORUM-SENSING ARCHITECTURE IN ROUXIELLA BADENSIS
# ============================================================
#
# Main question:
#
# Does Rouxiella badensis possess an ExpI/ExpR-like AHL
# quorum-sensing system comparable to the system regulating
# solanimycin in Dickeya solani?
#
#
# Reference:
#
#   Dickeya solani MK10
#   GCA_000365285.1
#
# MK10 is used because solanimycin and ExpIR regulation were
# experimentally investigated in this strain.
#
#
# Primary Rouxiella strain:
#
#   Rouxiella badensis 20GA0316
#   GCF_020740305.1
#
#
# Current workflow:
#
#   STEP 01 - Create directories / verify environment
#   STEP 02 - Download D. solani MK10
#   STEP 03 - Identify ExpI/LuxI and ExpR/LuxR candidates in MK10
#   STEP 04 - Create BLAST database for R. badensis 20GA0316
#   STEP 05 - Search MK10 ExpI candidates against 20GA0316
#   STEP 06 - Search MK10 ExpR candidates against 20GA0316
#   STEP 07 - Summarize best hits
#
#
# CHECKPOINT SYSTEM:
#
# Each completed step creates:
#
#   /scratch/al98750/Roux/05_QS_ExpIR/.state/STEP.done
#
# If this script is run again, completed steps are skipped.
#
# Therefore future steps can simply be added below this script.
#
# ============================================================



# ============================================================
# PROJECT PATHS
# ============================================================

ROOT="/scratch/al98750/Roux"

WORK="${ROOT}/05_QS_ExpIR"

REF="${WORK}/01_Dsolani_MK10"

QUERY="${WORK}/02_queries"

BLAST_DB="${WORK}/03_BLAST_database"

RAW_RESULTS="${WORK}/04_raw_results"

RESULTS="${WORK}/05_results"

CONTEXT="${WORK}/06_genomic_context"

STATE="${WORK}/.state"

TMP="${WORK}/tmp"

THREADS="${SLURM_CPUS_PER_TASK:-4}"


# ------------------------------------------------------------
# Primary R. badensis genome
# ------------------------------------------------------------

RB_ACC="GCF_020740305.1"

RB_GENOME="${ROOT}/02_genomes/${RB_ACC}.fna"


# ------------------------------------------------------------
# D. solani MK10 reference
# ------------------------------------------------------------

MK10_ACC="GCA_000365285.1"



# ============================================================
# CREATE ALL OUTPUT DIRECTORIES
# ============================================================

mkdir -p \
    "${WORK}" \
    "${REF}" \
    "${QUERY}" \
    "${BLAST_DB}" \
    "${RAW_RESULTS}" \
    "${RESULTS}" \
    "${CONTEXT}" \
    "${STATE}" \
    "${TMP}"



# ============================================================
# HELPER FUNCTIONS
# ============================================================

step_done () {

    [[ -f "${STATE}/$1.done" ]]

}


mark_done () {

    date --iso-8601=seconds > "${STATE}/$1.done"

}


die () {

    echo
    echo "ERROR: $*" >&2
    echo
    exit 1

}



# ============================================================
# JOB INFORMATION
# ============================================================

echo
echo "============================================================"
echo " ROUXIELLA ExpIR DISCOVERY PIPELINE"
echo "============================================================"
echo

echo "Job ID:            ${SLURM_JOB_ID:-interactive}"
echo "Node:              $(hostname)"
echo "Threads:           ${THREADS}"
echo

echo "R. badensis:"
echo "  ${RB_ACC}"
echo

echo "D. solani reference:"
echo "  ${MK10_ACC}"
echo

echo "Working directory:"
echo "  ${WORK}"
echo

echo "Started:"
echo "  $(date)"
echo



# ============================================================
# STEP 01
# ACTIVATE SOFTWARE ENVIRONMENT AND VERIFY INPUTS
# ============================================================

if ! step_done "01_environment"; then

    echo
    echo "============================================================"
    echo "STEP 01: Environment and input checks"
    echo "============================================================"
    echo


    source "$(conda info --base)/etc/profile.d/conda.sh"

    conda activate "${ROOT}/envs/roux"


    for PROGRAM in \
        datasets \
        seqkit \
        makeblastdb \
        tblastn \
        unzip

    do

        command -v "${PROGRAM}" >/dev/null 2>&1 || \
            die "${PROGRAM} was not found."

        echo "FOUND: ${PROGRAM}"

    done


    [[ -s "${RB_GENOME}" ]] || \
        die "20GA0316 genome not found: ${RB_GENOME}"


    echo
    echo "Rouxiella genome:"
    echo "${RB_GENOME}"
    echo


    mark_done "01_environment"

else

    echo
    echo ">>> STEP 01 already completed — skipping"

fi



# ============================================================
# ACTIVATE ENVIRONMENT FOR ALL SUBSEQUENT STEPS
# ============================================================

source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate "${ROOT}/envs/roux"



# ============================================================
# STEP 02
# DOWNLOAD DICKEYA SOLANI MK10
# ============================================================

if ! step_done "02_download_MK10"; then

    echo
    echo "============================================================"
    echo "STEP 02: Downloading Dickeya solani MK10"
    echo "============================================================"
    echo


    MK10_ZIP="${REF}/MK10_NCBI.zip"

    MK10_PACKAGE="${REF}/package"


    rm -f "${MK10_ZIP}"

    rm -rf "${MK10_PACKAGE}"


    datasets download genome accession \
        "${MK10_ACC}" \
        --include genome,gff3,gbff,protein,cds,rna \
        --filename "${MK10_ZIP}"


    [[ -s "${MK10_ZIP}" ]] || \
        die "MK10 NCBI download failed."


    mkdir -p "${MK10_PACKAGE}"


    unzip -q \
        "${MK10_ZIP}" \
        -d "${MK10_PACKAGE}"


    MK10_DIR="${MK10_PACKAGE}/ncbi_dataset/data/${MK10_ACC}"


    [[ -d "${MK10_DIR}" ]] || \
        die "MK10 NCBI package directory was not found."


    mark_done "02_download_MK10"

else

    echo
    echo ">>> STEP 02 already completed — skipping"

fi



# ============================================================
# DEFINE MK10 FILES
# ============================================================

MK10_DIR="${REF}/package/ncbi_dataset/data/${MK10_ACC}"


MK10_GFF="$(find "${MK10_DIR}" \
    -maxdepth 1 \
    -type f \
    \( -name "*.gff" -o -name "*.gff3" \) \
    | head -1)"


MK10_FAA="$(find "${MK10_DIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.faa" \
    | head -1)"


MK10_FNA="$(find "${MK10_DIR}" \
    -maxdepth 1 \
    -type f \
    -name "*genomic.fna" \
    | head -1)"


[[ -n "${MK10_GFF}" && -s "${MK10_GFF}" ]] || \
    die "MK10 GFF file was not found."


[[ -n "${MK10_FAA}" && -s "${MK10_FAA}" ]] || \
    die "MK10 protein FASTA was not found."


[[ -n "${MK10_FNA}" && -s "${MK10_FNA}" ]] || \
    die "MK10 genome FASTA was not found."



# ============================================================
# STEP 03
# IDENTIFY ExpI / ExpR CANDIDATES IN MK10
# ============================================================

if ! step_done "03_find_MK10_ExpIR"; then

    echo
    echo "============================================================"
    echo "STEP 03: Finding ExpI / ExpR candidates in MK10"
    echo "============================================================"
    echo


    # --------------------------------------------------------
    # Save broad annotation hits for manual inspection.
    # --------------------------------------------------------

    grep -iE \
'expI|expR|luxI|luxR|homoserine.lactone|AHL|autoinducer|quorum' \
        "${MK10_GFF}" \
        > "${RESULTS}/MK10_QS_annotation_hits.txt" \
        || true


    # --------------------------------------------------------
    # Search FASTA headers.
    # --------------------------------------------------------

    grep '^>' "${MK10_FAA}" \
        | grep -iE \
'expI|expR|luxI|luxR|homoserine.lactone|autoinducer|quorum' \
        > "${RESULTS}/MK10_QS_protein_headers.txt" \
        || true


    # --------------------------------------------------------
    # EXP I / LUX I CANDIDATES
    #
    # Broad patterns are intentional because annotation names
    # may vary.
    # --------------------------------------------------------

    seqkit grep \
        -r \
        -i \
        -p 'ExpI|LuxI|homoserine.lactone.*synthase|autoinducer.*synthase' \
        "${MK10_FAA}" \
        > "${QUERY}/MK10_ExpI_candidates.faa" \
        || true


    # --------------------------------------------------------
    # EXP R / LUX R CANDIDATES
    # --------------------------------------------------------

    seqkit grep \
        -r \
        -i \
        -p 'ExpR|LuxR|quorum.*regulator|autoinducer.*regulator' \
        "${MK10_FAA}" \
        > "${QUERY}/MK10_ExpR_candidates.faa" \
        || true


    echo
    echo "ExpI/LuxI candidate proteins:"
    echo

    grep '^>' \
        "${QUERY}/MK10_ExpI_candidates.faa" \
        || true


    echo
    echo "ExpR/LuxR candidate proteins:"
    echo

    grep '^>' \
        "${QUERY}/MK10_ExpR_candidates.faa" \
        || true


    echo


    mark_done "03_find_MK10_ExpIR"

else

    echo
    echo ">>> STEP 03 already completed — skipping"

fi



# ============================================================
# STEP 04
# CREATE BLAST DATABASE FOR R. BADENSIS 20GA0316
# ============================================================

if ! step_done "04_Rb_BLAST_database"; then

    echo
    echo "============================================================"
    echo "STEP 04: Building R. badensis BLAST database"
    echo "============================================================"
    echo


    rm -f "${BLAST_DB}/Rb20GA0316."*


    makeblastdb \
        -in "${RB_GENOME}" \
        -dbtype nucl \
        -parse_seqids \
        -out "${BLAST_DB}/Rb20GA0316"


    mark_done "04_Rb_BLAST_database"

else

    echo
    echo ">>> STEP 04 already completed — skipping"

fi



# ============================================================
# STEP 05
# SEARCH ExpI / LuxI AGAINST R. BADENSIS 20GA0316
# ============================================================

if ! step_done "05_ExpI_search"; then

    echo
    echo "============================================================"
    echo "STEP 05: Searching for ExpI/LuxI in 20GA0316"
    echo "============================================================"
    echo


    OUT="${RAW_RESULTS}/ExpI_vs_Rb20GA0316.tsv"


    if [[ -s "${QUERY}/MK10_ExpI_candidates.faa" ]]; then


        tblastn \
            -query "${QUERY}/MK10_ExpI_candidates.faa" \
            -db "${BLAST_DB}/Rb20GA0316" \
            -evalue 1e-10 \
            -max_target_seqs 50 \
            -num_threads "${THREADS}" \
            -outfmt \
'6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore' \
            > "${OUT}"


    else


        echo \
"WARNING: No MK10 ExpI candidate was automatically extracted." \
            >&2


        : > "${OUT}"


    fi


    mark_done "05_ExpI_search"

else

    echo
    echo ">>> STEP 05 already completed — skipping"

fi



# ============================================================
# STEP 06
# SEARCH ExpR / LuxR AGAINST R. BADENSIS 20GA0316
# ============================================================

if ! step_done "06_ExpR_search"; then

    echo
    echo "============================================================"
    echo "STEP 06: Searching for ExpR/LuxR in 20GA0316"
    echo "============================================================"
    echo


    OUT="${RAW_RESULTS}/ExpR_vs_Rb20GA0316.tsv"


    if [[ -s "${QUERY}/MK10_ExpR_candidates.faa" ]]; then


        tblastn \
            -query "${QUERY}/MK10_ExpR_candidates.faa" \
            -db "${BLAST_DB}/Rb20GA0316" \
            -evalue 1e-10 \
            -max_target_seqs 100 \
            -num_threads "${THREADS}" \
            -outfmt \
'6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore' \
            > "${OUT}"


    else


        echo \
"WARNING: No MK10 ExpR candidate was automatically extracted." \
            >&2


        : > "${OUT}"


    fi


    mark_done "06_ExpR_search"

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi



# ============================================================
# STEP 07
# CREATE READABLE HIT SUMMARIES
# ============================================================

if ! step_done "07_summary"; then

    echo
    echo "============================================================"
    echo "STEP 07: Creating ExpI / ExpR summary tables"
    echo "============================================================"
    echo


    HEADER=$'query\ttarget\tidentity_pct\talignment_aa\tquery_length\tmismatches\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'


    # --------------------------------------------------------
    # ExpI
    # --------------------------------------------------------

    echo "${HEADER}" \
        > "${RESULTS}/ExpI_best_hits.tsv"


    if [[ -s "${RAW_RESULTS}/ExpI_vs_Rb20GA0316.tsv" ]]; then

        sort \
            -t $'\t' \
            -k13,13nr \
            "${RAW_RESULTS}/ExpI_vs_Rb20GA0316.tsv" \
            | head -20 \
            >> "${RESULTS}/ExpI_best_hits.tsv"

    fi


    # --------------------------------------------------------
    # ExpR
    # --------------------------------------------------------

    echo "${HEADER}" \
        > "${RESULTS}/ExpR_best_hits.tsv"


    if [[ -s "${RAW_RESULTS}/ExpR_vs_Rb20GA0316.tsv" ]]; then

        sort \
            -t $'\t' \
            -k13,13nr \
            "${RAW_RESULTS}/ExpR_vs_Rb20GA0316.tsv" \
            | head -40 \
            >> "${RESULTS}/ExpR_best_hits.tsv"

    fi


    # --------------------------------------------------------
    # Simple candidate counts
    # --------------------------------------------------------

    {

        echo -e \
"analysis\tnumber_of_BLAST_HSPs"


        EXP_I_COUNT=$(grep -vc '^$' \
            "${RAW_RESULTS}/ExpI_vs_Rb20GA0316.tsv" \
            || true)


        EXP_R_COUNT=$(grep -vc '^$' \
            "${RAW_RESULTS}/ExpR_vs_Rb20GA0316.tsv" \
            || true)


        echo -e \
"ExpI_vs_20GA0316\t${EXP_I_COUNT}"


        echo -e \
"ExpR_vs_20GA0316\t${EXP_R_COUNT}"


    } > "${RESULTS}/ExpIR_hit_counts.tsv"


    mark_done "07_summary"

else

    echo
    echo ">>> STEP 07 already completed — skipping"

fi



# ============================================================
# FINAL REPORT
# ============================================================

echo
echo
echo "============================================================"
echo " ExpIR DISCOVERY PIPELINE COMPLETE"
echo "============================================================"
echo

echo "Working directory:"
echo "  ${WORK}"
echo

echo "MK10 QS annotation hits:"
echo "  ${RESULTS}/MK10_QS_annotation_hits.txt"
echo

echo "MK10 QS protein headers:"
echo "  ${RESULTS}/MK10_QS_protein_headers.txt"
echo

echo "MK10 ExpI candidate proteins:"
echo "  ${QUERY}/MK10_ExpI_candidates.faa"
echo

echo "MK10 ExpR candidate proteins:"
echo "  ${QUERY}/MK10_ExpR_candidates.faa"
echo

echo "ExpI hits in 20GA0316:"
echo "  ${RESULTS}/ExpI_best_hits.tsv"
echo

echo "ExpR hits in 20GA0316:"
echo "  ${RESULTS}/ExpR_best_hits.tsv"
echo

echo "Hit counts:"
echo "  ${RESULTS}/ExpIR_hit_counts.tsv"
echo

echo "Completed checkpoints:"
echo "  ${STATE}/"
echo

echo "Finished:"
echo "  $(date)"
echo



# ============================================================
# ADD FUTURE STEPS BELOW THIS POINT
# ============================================================
#
# IMPORTANT:
#
# Do NOT modify the names of completed checkpoint steps above.
#
# For example, later we can append:
#
#
# if ! step_done "08_extract_Rb_candidates"; then
#
#     echo "STEP 08: Extract exact R. badensis candidate genes"
#
#     COMMANDS
#
#     mark_done "08_extract_Rb_candidates"
#
# else
#
#     echo "STEP 08 already completed — skipping"
#
# fi
#
#
# Then:
#
#   STEP 09 - identify genes in Bakta annotation
#   STEP 10 - extract +/- 20 kb genomic neighborhood
#   STEP 11 - locate sol cluster
#   STEP 12 - calculate ExpIR-to-sol genomic distance
#   STEP 13 - compare all R. badensis genomes
#   STEP 14 - compare architecture with D. solani MK10
#
# Re-running this script will skip STEP 01-07 because their
# .done files already exist.
#
# ============================================================