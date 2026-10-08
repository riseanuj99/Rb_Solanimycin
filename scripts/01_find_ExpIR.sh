#!/usr/bin/env bash

#SBATCH --job-name=Roux_ExpIR
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/Roux_ExpIR_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/Roux_ExpIR_%j.err

set -euo pipefail
shopt -s nullglob


# ============================================================
# OBJECTIVE 1
# AHL QUORUM-SENSING ARCHITECTURE IN ROUXIELLA BADENSIS
# ============================================================
#
# Question:
#
# Does Rouxiella badensis possess an ExpI/ExpR-like AHL
# quorum-sensing system comparable to that regulating
# solanimycin in Dickeya solani?
#
#
# Reference:
#
#   Dickeya solani MK10
#   GCA_000365285.1
#
#
# Primary target:
#
#   Rouxiella badensis 20GA0316
#   GCF_020740305.1
#
#
# Workflow:
#
#   01 - Check environment
#   02 - Download D. solani MK10 genome
#   03 - Annotate MK10 consistently with Bakta
#   04 - Identify ExpI/LuxI and ExpR/LuxR candidates in MK10
#   05 - Build 20GA0316 BLAST database
#   06 - Search ExpI candidates against 20GA0316
#   07 - Search ExpR candidates against 20GA0316
#   08 - Create summary tables
#
#
# Checkpoints:
#
#   /scratch/al98750/Roux/05_QS_ExpIR/.state/
#
# Completed steps are NOT rerun.
#
# Future STEP 09, STEP 10, etc. can be appended later.
#
# ============================================================



# ============================================================
# PATHS
# ============================================================

ROOT="/scratch/al98750/Roux"

WORK="${ROOT}/05_QS_ExpIR"

REF="${WORK}/01_Dsolani_MK10"

MK10_ANNOT="${REF}/bakta_annotation"

QUERY="${WORK}/02_queries"

BLAST_DB="${WORK}/03_BLAST_database"

RAW_RESULTS="${WORK}/04_raw_results"

RESULTS="${WORK}/05_results"

CONTEXT="${WORK}/06_genomic_context"

STATE="${WORK}/.state"

TMP="${WORK}/tmp"


THREADS="${SLURM_CPUS_PER_TASK:-4}"


# ============================================================
# ACCESSIONS
# ============================================================

RB_ACC="GCF_020740305.1"

RB_GENOME="${ROOT}/02_genomes/${RB_ACC}.fna"


MK10_ACC="GCA_000365285.1"


# ============================================================
# BAKTA DATABASE
# ============================================================

BAKTA_DB="${ROOT}/db/bakta/db-light"



# ============================================================
# CREATE ALL NEEDED DIRECTORIES
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

echo "Job ID:"
echo "  ${SLURM_JOB_ID:-interactive}"
echo

echo "Node:"
echo "  $(hostname)"
echo

echo "Threads:"
echo "  ${THREADS}"
echo

echo "R. badensis:"
echo "  ${RB_ACC}"
echo

echo "D. solani reference:"
echo "  ${MK10_ACC}"
echo

echo "Work directory:"
echo "  ${WORK}"
echo

echo "Started:"
echo "  $(date)"
echo



# ============================================================
# ACTIVATE EXISTING ROUX ENVIRONMENT
# ============================================================

source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate "${ROOT}/envs/roux"



# ============================================================
# STEP 01
# CHECK ENVIRONMENT
# ============================================================

if ! step_done "01_environment"; then

    echo
    echo "============================================================"
    echo "STEP 01: Environment and input checks"
    echo "============================================================"
    echo


    for PROGRAM in \
        datasets \
        bakta \
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
        die "20GA0316 genome missing: ${RB_GENOME}"


    [[ -d "${BAKTA_DB}" ]] || \
        die "Bakta database missing: ${BAKTA_DB}"


    mark_done "01_environment"

else

    echo
    echo ">>> STEP 01 already completed — skipping"

fi



# ============================================================
# ALWAYS VERIFY REQUIRED SOFTWARE
# ============================================================

for PROGRAM in \
    datasets \
    bakta \
    seqkit \
    makeblastdb \
    tblastn \
    unzip

do

    command -v "${PROGRAM}" >/dev/null 2>&1 || \
        die "${PROGRAM} was not found."

done



# ============================================================
# STEP 02
# DOWNLOAD D. SOLANI MK10
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


    # --------------------------------------------------------
    # Genome sequence is the critical requirement.
    #
    # Annotation files may not be provided for this GenBank
    # assembly, so MK10 will be annotated locally with Bakta
    # in STEP 03.
    # --------------------------------------------------------

    datasets download genome accession \
        "${MK10_ACC}" \
        --include genome \
        --filename "${MK10_ZIP}"


    [[ -s "${MK10_ZIP}" ]] || \
        die "MK10 download failed."


    mkdir -p "${MK10_PACKAGE}"


    unzip -q \
        "${MK10_ZIP}" \
        -d "${MK10_PACKAGE}"


    mark_done "02_download_MK10"

else

    echo
    echo ">>> STEP 02 already completed — skipping"

fi



# ============================================================
# LOCATE MK10 GENOME
# ============================================================

MK10_PACKAGE="${REF}/package"


MK10_FNA="$(find "${MK10_PACKAGE}/ncbi_dataset/data" \
    -type f \
    -name "*.fna" \
    | head -1)"


[[ -n "${MK10_FNA}" && -s "${MK10_FNA}" ]] || {

    echo
    echo "Files currently present in MK10 package:"
    find "${MK10_PACKAGE}" -type f | sort
    echo

    die "MK10 genome FASTA was not found."

}


echo
echo "MK10 genome:"
echo "  ${MK10_FNA}"
echo



# ============================================================
# STEP 03
# ANNOTATE MK10 WITH BAKTA
# ============================================================
#
# This solves the previous failure.
#
# GCA_000365285.1 does not provide the annotation files that
# our previous script expected.
#
# We therefore annotate MK10 ourselves using exactly the same
# Bakta installation/database used for R. badensis.
#
# This also makes comparisons more consistent.
#
# ============================================================

if ! step_done "03_annotate_MK10"; then

    echo
    echo "============================================================"
    echo "STEP 03: Annotating D. solani MK10 with Bakta"
    echo "============================================================"
    echo


    if [[ -d "${MK10_ANNOT}" ]]; then

        echo "Removing incomplete previous MK10 annotation:"
        echo "  ${MK10_ANNOT}"

        rm -rf "${MK10_ANNOT}"

    fi


    # IMPORTANT:
    # Do not mkdir MK10_ANNOT.
    # Bakta creates the output directory itself.


    bakta \
        --db "${BAKTA_DB}" \
        --output "${MK10_ANNOT}" \
        --prefix "Dsolani_MK10" \
        --genus Dickeya \
        --species solani \
        --gram - \
        --threads "${THREADS}" \
        --keep-contig-headers \
        "${MK10_FNA}"


    [[ -s "${MK10_ANNOT}/Dsolani_MK10.gff3" ]] || \
        die "MK10 Bakta GFF3 was not created."


    [[ -s "${MK10_ANNOT}/Dsolani_MK10.faa" ]] || \
        die "MK10 Bakta protein FASTA was not created."


    [[ -s "${MK10_ANNOT}/Dsolani_MK10.tsv" ]] || \
        die "MK10 Bakta TSV was not created."


    mark_done "03_annotate_MK10"

else

    echo
    echo ">>> STEP 03 already completed — skipping"

fi



# ============================================================
# DEFINE STANDARDIZED MK10 ANNOTATION FILES
# ============================================================

MK10_GFF="${MK10_ANNOT}/Dsolani_MK10.gff3"

MK10_FAA="${MK10_ANNOT}/Dsolani_MK10.faa"

MK10_TSV="${MK10_ANNOT}/Dsolani_MK10.tsv"

MK10_FFN="${MK10_ANNOT}/Dsolani_MK10.ffn"


[[ -s "${MK10_GFF}" ]] || \
    die "MK10 GFF missing: ${MK10_GFF}"


[[ -s "${MK10_FAA}" ]] || \
    die "MK10 FAA missing: ${MK10_FAA}"


[[ -s "${MK10_TSV}" ]] || \
    die "MK10 TSV missing: ${MK10_TSV}"



# ============================================================
# STEP 04
# IDENTIFY ExpI / LuxI AND ExpR / LuxR CANDIDATES IN MK10
# ============================================================

if ! step_done "04_find_MK10_ExpIR"; then

    echo
    echo "============================================================"
    echo "STEP 04: Finding ExpI/LuxI and ExpR/LuxR in MK10"
    echo "============================================================"
    echo


    # --------------------------------------------------------
    # Broad TSV search
    # --------------------------------------------------------

    grep -iE \
'expI|expR|luxI|luxR|homoserine lactone|autoinducer|quorum' \
        "${MK10_TSV}" \
        > "${RESULTS}/MK10_QS_annotation_hits.txt" \
        || true


    # --------------------------------------------------------
    # Protein FASTA headers
    # --------------------------------------------------------

    grep '^>' "${MK10_FAA}" \
        | grep -iE \
'expI|expR|luxI|luxR|homoserine.lactone|autoinducer|quorum' \
        > "${RESULTS}/MK10_QS_protein_headers.txt" \
        || true


    # ========================================================
    # ExpI / LuxI candidate proteins
    # ========================================================

    seqkit grep \
        -r \
        -i \
        -p \
'ExpI|LuxI|acyl.*homoserine.*lactone.*synthase|homoserine.*lactone.*synthase|autoinducer.*synthase' \
        "${MK10_FAA}" \
        > "${QUERY}/MK10_ExpI_candidates.faa" \
        || true


    # ========================================================
    # ExpR / LuxR candidate proteins
    # ========================================================

    seqkit grep \
        -r \
        -i \
        -p \
'ExpR|LuxR|LuxR.*family|quorum.*regulator|autoinducer.*regulator' \
        "${MK10_FAA}" \
        > "${QUERY}/MK10_ExpR_candidates.faa" \
        || true


    echo
    echo "MK10 ExpI/LuxI candidates:"
    echo

    grep '^>' \
        "${QUERY}/MK10_ExpI_candidates.faa" \
        || true


    echo
    echo "MK10 ExpR/LuxR candidates:"
    echo

    grep '^>' \
        "${QUERY}/MK10_ExpR_candidates.faa" \
        || true


    echo


    mark_done "04_find_MK10_ExpIR"

else

    echo
    echo ">>> STEP 04 already completed — skipping"

fi



# ============================================================
# STEP 05
# CREATE 20GA0316 BLAST DATABASE
# ============================================================

if ! step_done "05_Rb_BLAST_database"; then

    echo
    echo "============================================================"
    echo "STEP 05: Building 20GA0316 BLAST database"
    echo "============================================================"
    echo


    rm -f "${BLAST_DB}/Rb20GA0316."*


    makeblastdb \
        -in "${RB_GENOME}" \
        -dbtype nucl \
        -parse_seqids \
        -out "${BLAST_DB}/Rb20GA0316"


    mark_done "05_Rb_BLAST_database"

else

    echo
    echo ">>> STEP 05 already completed — skipping"

fi



# ============================================================
# STEP 06
# SEARCH ExpI / LuxI AGAINST 20GA0316
# ============================================================

if ! step_done "06_ExpI_search"; then

    echo
    echo "============================================================"
    echo "STEP 06: Searching ExpI/LuxI against 20GA0316"
    echo "============================================================"
    echo


    OUT="${RAW_RESULTS}/ExpI_vs_Rb20GA0316.tsv"


    : > "${OUT}"


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
"WARNING: No ExpI candidate automatically identified in MK10." \
            >&2

    fi


    mark_done "06_ExpI_search"

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi



# ============================================================
# STEP 07
# SEARCH ExpR / LuxR AGAINST 20GA0316
# ============================================================

if ! step_done "07_ExpR_search"; then

    echo
    echo "============================================================"
    echo "STEP 07: Searching ExpR/LuxR against 20GA0316"
    echo "============================================================"
    echo


    OUT="${RAW_RESULTS}/ExpR_vs_Rb20GA0316.tsv"


    : > "${OUT}"


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
"WARNING: No ExpR candidate automatically identified in MK10." \
            >&2

    fi


    mark_done "07_ExpR_search"

else

    echo
    echo ">>> STEP 07 already completed — skipping"

fi



# ============================================================
# STEP 08
# SUMMARIZE RESULTS
# ============================================================

if ! step_done "08_summary"; then

    echo
    echo "============================================================"
    echo "STEP 08: Creating summary tables"
    echo "============================================================"
    echo


    HEADER=$'query\ttarget\tidentity_pct\talignment_aa\tquery_length\tmismatches\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore'


    # ========================================================
    # ExpI
    # ========================================================

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


    # ========================================================
    # ExpR
    # ========================================================

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


    # ========================================================
    # Counts
    # ========================================================

    EXP_I_COUNT=$(

        awk 'NF > 0 {n++} END {print n+0}' \
            "${RAW_RESULTS}/ExpI_vs_Rb20GA0316.tsv"

    )


    EXP_R_COUNT=$(

        awk 'NF > 0 {n++} END {print n+0}' \
            "${RAW_RESULTS}/ExpR_vs_Rb20GA0316.tsv"

    )


    {

        echo -e \
"analysis\tBLAST_HSP_count"


        echo -e \
"ExpI_vs_20GA0316\t${EXP_I_COUNT}"


        echo -e \
"ExpR_vs_20GA0316\t${EXP_R_COUNT}"


    } > "${RESULTS}/ExpIR_hit_counts.tsv"


    mark_done "08_summary"

else

    echo
    echo ">>> STEP 08 already completed — skipping"

fi



# ============================================================
# FINAL REPORT
# ============================================================

echo
echo
echo "============================================================"
echo " ExpIR DISCOVERY COMPLETE"
echo "============================================================"
echo

echo "MK10 Bakta annotation:"
echo "  ${MK10_ANNOT}"
echo

echo "MK10 QS annotation hits:"
echo "  ${RESULTS}/MK10_QS_annotation_hits.txt"
echo

echo "MK10 ExpI candidate FASTA:"
echo "  ${QUERY}/MK10_ExpI_candidates.faa"
echo

echo "MK10 ExpR candidate FASTA:"
echo "  ${QUERY}/MK10_ExpR_candidates.faa"
echo

echo "ExpI best hits:"
echo "  ${RESULTS}/ExpI_best_hits.tsv"
echo

echo "ExpR best hits:"
echo "  ${RESULTS}/ExpR_best_hits.tsv"
echo

echo "Hit counts:"
echo "  ${RESULTS}/ExpIR_hit_counts.tsv"
echo

echo "Finished:"
echo "  $(date)"
echo



# ============================================================
# FUTURE STEPS
# ============================================================
#
# Add future analysis here without changing checkpoint names:
#
#
# STEP 09:
#   map Rouxiella BLAST hits to Bakta genes
#
# STEP 10:
#   identify exact ExpI / ExpR orthologs
#
# STEP 11:
#   extract +/- 20 kb neighborhoods
#
# STEP 12:
#   identify sol cluster in 20GA0316
#
# STEP 13:
#   calculate ExpIR-to-sol distance
#
# STEP 14:
#   compare all 17 R. badensis genomes
#
# STEP 15:
#   compare architecture to D. solani MK10
#
# ============================================================