#!/usr/bin/env bash

#SBATCH --job-name=Roux_annot
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/Roux_annot_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/Roux_annot_%j.err

set -euo pipefail
shopt -s nullglob


# ============================================================
# ROUXIELLA BADENSIS GENOME + ANNOTATION PIPELINE
# ============================================================
#
# Main project:
#   /scratch/al98750/Roux
#
# Current goals:
#
#   00. Create software environment
#   01. Query NCBI for all current Rouxiella badensis genomes
#   02. Download genomes + existing NCBI annotations
#   03. Organize genome files
#   04. Download Bakta database
#   05. Reannotate every genome with Bakta
#   06. Summarize CDS / tRNA / rRNA / ncRNA / sRNA features
#   07. Search specifically for ArcZ using the Rfam model
#
# Checkpoint design:
#
#   Completed steps get:
#       /scratch/al98750/Roux/.state/STEP.done
#
#   Individual genome annotations also get their own .done files.
#
# Therefore:
#   re-running this script DOES NOT repeat completed work.
#
# Future analyses can simply be added as STEP 08, STEP 09, etc.
#
# ============================================================



# ============================================================
# PROJECT SETTINGS
# ============================================================

ROOT="/scratch/al98750/Roux"

THREADS="${SLURM_CPUS_PER_TASK:-8}"


# ------------------------------------------------------------
# Main directories
# ------------------------------------------------------------

SCRIPTS="${ROOT}/scripts"

LOGS="${ROOT}/logs"

STATE="${ROOT}/.state"

ENV_ROOT="${ROOT}/envs"

ENV_DIR="${ENV_ROOT}/roux"


# ------------------------------------------------------------
# NCBI
# ------------------------------------------------------------

NCBI="${ROOT}/01_NCBI"

NCBI_PACKAGE="${NCBI}/package"

NCBI_ORIGINAL="${NCBI}/original_annotations"


# ------------------------------------------------------------
# Genome FASTA files
# ------------------------------------------------------------

GENOMES="${ROOT}/02_genomes"


# ------------------------------------------------------------
# Annotation
# ------------------------------------------------------------

ANNOT_ROOT="${ROOT}/03_annotation"

BAKTA_ANNOT="${ANNOT_ROOT}/bakta"


# ------------------------------------------------------------
# Results
# ------------------------------------------------------------

RESULTS="${ROOT}/04_results"

GENOME_QC="${RESULTS}/genome_QC"

FEATURE_RESULTS="${RESULTS}/annotation_features"

ARCS="${RESULTS}/ArcZ"


# ------------------------------------------------------------
# Databases
# ------------------------------------------------------------

DB_ROOT="${ROOT}/db"

BAKTA_DB_ROOT="${DB_ROOT}/bakta"


# ============================================================
# BAKTA DATABASE TYPE
# ============================================================
#
# "light" is recommended initially.
#
# IMPORTANT:
# Bakta LIGHT still contains the complete non-coding feature
# database used for:
#
#   ncRNA / sRNA
#   tRNA
#   tmRNA
#   rRNA
#   riboswitches
#   RNA leaders
#   CRISPR
#   etc.
#
# CDS/ORF PREDICTION IS STILL PERFORMED.
#
# FULL mainly provides more detailed protein functional
# annotation but requires ~84 GB when unpacked.
#
# ============================================================

BAKTA_DB_TYPE="light"


if [[ "${BAKTA_DB_TYPE}" == "light" ]]; then

    BAKTA_DB="${BAKTA_DB_ROOT}/db-light"

else

    BAKTA_DB="${BAKTA_DB_ROOT}/db"

fi



# ============================================================
# CREATE ALL PROJECT DIRECTORIES
# ============================================================

mkdir -p \
    "${SCRIPTS}" \
    "${LOGS}" \
    "${STATE}" \
    "${ENV_ROOT}" \
    "${NCBI}" \
    "${NCBI_PACKAGE}" \
    "${NCBI_ORIGINAL}" \
    "${GENOMES}" \
    "${ANNOT_ROOT}" \
    "${BAKTA_ANNOT}" \
    "${RESULTS}" \
    "${GENOME_QC}" \
    "${FEATURE_RESULTS}" \
    "${ARCS}" \
    "${DB_ROOT}" \
    "${BAKTA_DB_ROOT}"



# ============================================================
# LOG BASIC JOB INFORMATION
# ============================================================

echo
echo "============================================================"
echo " ROUXIELLA BADENSIS PIPELINE"
echo "============================================================"
echo
echo "Project:       ${ROOT}"
echo "Job ID:        ${SLURM_JOB_ID:-interactive}"
echo "Node:          $(hostname)"
echo "Threads:       ${THREADS}"
echo "Start time:    $(date)"
echo
echo "============================================================"
echo



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
# STEP 00
# CREATE SOFTWARE ENVIRONMENT
# ============================================================
#
# Software:
#
#   NCBI Datasets
#   Bakta
#   BLAST+
#   Infernal
#   SeqKit
#   jq
#   wget
#   unzip
#
# ============================================================

if ! step_done "00_environment"; then

    echo
    echo ">>> STEP 00: Setting up software environment"
    echo


    if ! command -v conda >/dev/null 2>&1; then

        die "conda is not available in this Slurm job."

    fi


    CONDA_BASE="$(conda info --base)"

    source "${CONDA_BASE}/etc/profile.d/conda.sh"


    if command -v mamba >/dev/null 2>&1; then

        SOLVER="mamba"

    else

        SOLVER="conda"

    fi


    PACKAGES=(

        ncbi-datasets-cli
        bakta
        blast
        infernal
        seqkit
        jq
        wget
        unzip

    )


    if [[ -d "${ENV_DIR}/conda-meta" ]]; then

        echo "Existing environment detected."
        echo "Checking/installing required packages."

        "${SOLVER}" install -y \
            -p "${ENV_DIR}" \
            -c conda-forge \
            -c bioconda \
            "${PACKAGES[@]}"

    else

        echo "Creating new Rouxiella environment."

        "${SOLVER}" create -y \
            -p "${ENV_DIR}" \
            -c conda-forge \
            -c bioconda \
            "${PACKAGES[@]}"

    fi


    mark_done "00_environment"

else

    echo
    echo ">>> STEP 00 already completed — skipping"

fi



# ============================================================
# ACTIVATE ENVIRONMENT
# ============================================================

CONDA_BASE="$(conda info --base)"

source "${CONDA_BASE}/etc/profile.d/conda.sh"

conda activate "${ENV_DIR}"



# ============================================================
# VERIFY IMPORTANT PROGRAMS
# ============================================================

echo
echo "Software check:"
echo


for PROGRAM in \
    datasets \
    dataformat \
    bakta \
    bakta_db \
    blastn \
    cmsearch \
    seqkit \
    jq \
    wget \
    unzip

do

    if command -v "${PROGRAM}" >/dev/null 2>&1; then

        echo "FOUND: ${PROGRAM}"

    else

        die "${PROGRAM} was not found."

    fi

done


echo
echo "Versions:"
echo

datasets version || true

bakta --version || true

blastn -version | head -1 || true

cmsearch -h | head -2 || true



# ============================================================
# STEP 01
# QUERY NCBI FOR ALL CURRENT ROUXIELLA BADENSIS GENOMES
# ============================================================
#
# Exact taxon:
#
#   Rouxiella badensis
#
# This produces:
#
#   metadata JSONL
#   metadata TSV
#   accession list
#   separate table showing any 20GA0316 entry
#
# ============================================================

if ! step_done "01_ncbi_query"; then

    echo
    echo "============================================================"
    echo "STEP 01: Querying NCBI for Rouxiella badensis"
    echo "============================================================"
    echo


    datasets summary genome taxon \
        "Rouxiella badensis" \
        --tax-exact-match \
        --as-json-lines \
        > "${NCBI}/R_badensis_metadata.jsonl"


    dataformat tsv genome \
        --inputfile "${NCBI}/R_badensis_metadata.jsonl" \
        --fields \
accession,organism-name,organism-infraspecific-strain,assminfo-name,assminfo-level,assminfo-release-date,annotinfo-status,source_database \
        > "${NCBI}/R_badensis_assemblies.tsv"


    awk -F '\t' '

        NR > 1 && $1 ~ /^GC[AF]_/ {
            print $1
        }

    ' "${NCBI}/R_badensis_assemblies.tsv" \
        | sort -u \
        > "${NCBI}/accessions.txt"


    N_ACCESSIONS=$(wc -l < "${NCBI}/accessions.txt")


    if [[ "${N_ACCESSIONS}" -eq 0 ]]; then

        die "NCBI query returned zero Rouxiella badensis assemblies."

    fi


    # --------------------------------------------------------
    # Make a small table specifically for strain 20GA0316
    # --------------------------------------------------------

    {

        head -1 "${NCBI}/R_badensis_assemblies.tsv"

        grep -F "20GA0316" \
            "${NCBI}/R_badensis_assemblies.tsv" \
            || true

    } > "${NCBI}/R_badensis_20GA0316.tsv"


    echo
    echo "NCBI genomes found: ${N_ACCESSIONS}"
    echo

    cat "${NCBI}/R_badensis_assemblies.tsv"

    echo


    mark_done "01_ncbi_query"

else

    echo
    echo ">>> STEP 01 already completed — skipping"

fi



# ============================================================
# STEP 02
# DOWNLOAD GENOMES + ORIGINAL NCBI ANNOTATIONS
# ============================================================
#
# Retains:
#
#   genome FASTA
#   protein FASTA
#   CDS nucleotide FASTA
#   RNA FASTA
#   GFF3
#   GenBank flatfile
#   sequence report
#
# ============================================================

if ! step_done "02_ncbi_download"; then

    echo
    echo "============================================================"
    echo "STEP 02: Downloading NCBI genomes + annotations"
    echo "============================================================"
    echo


    ZIP="${NCBI}/R_badensis_NCBI.zip"


    # Remove incomplete files from a previously interrupted run.

    rm -f "${ZIP}"

    rm -rf "${NCBI_PACKAGE}"

    mkdir -p "${NCBI_PACKAGE}"


    datasets download genome accession \
        --inputfile "${NCBI}/accessions.txt" \
        --include genome,gff3,gbff,rna,cds,protein,seq-report \
        --filename "${ZIP}"


    [[ -s "${ZIP}" ]] || \
        die "NCBI genome download did not create ${ZIP}"


    unzip -q \
        "${ZIP}" \
        -d "${NCBI_PACKAGE}"


    [[ -d "${NCBI_PACKAGE}/ncbi_dataset/data" ]] || \
        die "NCBI package was not unpacked correctly."


    mark_done "02_ncbi_download"

else

    echo
    echo ">>> STEP 02 already completed — skipping"

fi



# ============================================================
# STEP 03
# ORGANIZE GENOME FILES
# ============================================================

if ! step_done "03_prepare_genomes"; then

    echo
    echo "============================================================"
    echo "STEP 03: Organizing genome FASTA and NCBI annotation files"
    echo "============================================================"
    echo


    while read -r ACC; do


        [[ -n "${ACC}" ]] || continue


        SOURCE_DIR="${NCBI_PACKAGE}/ncbi_dataset/data/${ACC}"


        if [[ ! -d "${SOURCE_DIR}" ]]; then

            echo "WARNING: NCBI directory missing for ${ACC}"

            continue

        fi


        # ----------------------------------------------------
        # Find genomic FASTA
        # ----------------------------------------------------

        FNA=""


        CANDIDATES=(

            "${SOURCE_DIR}"/*_genomic.fna
            "${SOURCE_DIR}"/genomic.fna

        )


        if [[ "${#CANDIDATES[@]}" -gt 0 ]]; then

            FNA="${CANDIDATES[0]}"

        fi


        if [[ -z "${FNA}" || ! -s "${FNA}" ]]; then

            echo "WARNING: genome FASTA missing for ${ACC}"

            continue

        fi


        cp -f \
            "${FNA}" \
            "${GENOMES}/${ACC}.fna"


        # ----------------------------------------------------
        # Preserve original NCBI annotations
        # ----------------------------------------------------

        ORIGINAL_OUT="${NCBI_ORIGINAL}/${ACC}"

        mkdir -p "${ORIGINAL_OUT}"


        for FILE in \
            "${SOURCE_DIR}"/*.gff \
            "${SOURCE_DIR}"/*.gff3 \
            "${SOURCE_DIR}"/*.gbff \
            "${SOURCE_DIR}"/*.faa \
            "${SOURCE_DIR}"/*cds*.fna \
            "${SOURCE_DIR}"/*rna*.fna \
            "${SOURCE_DIR}"/*sequence_report*.jsonl

        do

            [[ -e "${FILE}" ]] || continue

            cp -f \
                "${FILE}" \
                "${ORIGINAL_OUT}/"

        done


        echo "Prepared ${ACC}"


    done < "${NCBI}/accessions.txt"


    GENOME_FILES=( "${GENOMES}"/*.fna )


    if [[ "${#GENOME_FILES[@]}" -eq 0 ]]; then

        die "No genome FASTA files were prepared."

    fi


    # --------------------------------------------------------
    # Genome statistics
    # --------------------------------------------------------

    seqkit stats \
        -T \
        "${GENOME_FILES[@]}" \
        > "${GENOME_QC}/genome_stats.tsv"


    echo
    echo "Genome FASTA files:"
    echo

    ls -lh "${GENOMES}"


    echo
    echo "Genome statistics:"
    echo

    column -t \
        "${GENOME_QC}/genome_stats.tsv" \
        || cat "${GENOME_QC}/genome_stats.tsv"


    mark_done "03_prepare_genomes"

else

    echo
    echo ">>> STEP 03 already completed — skipping"

fi



# ============================================================
# STEP 04
# DOWNLOAD BAKTA DATABASE
# ============================================================

if ! step_done "04_bakta_database" || \
   [[ ! -d "${BAKTA_DB}" ]]; then

    echo
    echo "============================================================"
    echo "STEP 04: Preparing Bakta ${BAKTA_DB_TYPE} database"
    echo "============================================================"
    echo


    if [[ "${BAKTA_DB_TYPE}" == "light" ]]; then

        EXPECTED_DB="${BAKTA_DB_ROOT}/db-light"

    else

        EXPECTED_DB="${BAKTA_DB_ROOT}/db"

    fi


    # If the step was interrupted before completion,
    # remove the incomplete target database.

    if ! step_done "04_bakta_database"; then

        rm -rf "${EXPECTED_DB}"

    fi


    if [[ ! -d "${EXPECTED_DB}" ]]; then

        bakta_db download \
            --output "${BAKTA_DB_ROOT}" \
            --type "${BAKTA_DB_TYPE}"

    fi


    [[ -d "${BAKTA_DB}" ]] || \
        die "Bakta database directory was not created."


    mark_done "04_bakta_database"

else

    echo
    echo ">>> STEP 04 already completed — skipping"

fi



# ============================================================
# STEP 05
# COMPREHENSIVE BAKTA ANNOTATION
# ============================================================
#
# Bakta identifies/predicts, among other features:
#
#   CDS
#   small ORFs
#   pseudogenes
#   tRNA
#   tmRNA
#   rRNA
#   ncRNA / known sRNA families
#   ncRNA regulatory regions
#   riboswitches
#   RNA leaders
#   CRISPR regions
#   origins
#
# Each genome gets an independent .done marker.
#
# If the job stops halfway through, the next run starts with
# the first unfinished genome.
#
# ============================================================

echo
echo "============================================================"
echo "STEP 05: Bakta annotation"
echo "============================================================"
echo


GENOME_FILES=( "${GENOMES}"/*.fna )


if [[ "${#GENOME_FILES[@]}" -eq 0 ]]; then

    die "No genome FASTA files are available for Bakta."

fi


for FNA in "${GENOME_FILES[@]}"; do


    ACC="$(basename "${FNA}" .fna)"


    OUT="${BAKTA_ANNOT}/${ACC}"


    GENOME_DONE="${OUT}/.annotation.done"


    if [[ -f "${GENOME_DONE}" ]]; then

        echo
        echo "${ACC}: already annotated — skipping"

        continue

    fi


    echo
    echo "------------------------------------------------------------"
    echo "ANNOTATING: ${ACC}"
    echo "------------------------------------------------------------"
    echo


    if [[ -d "${OUT}" ]]; then

    echo "Removing incomplete Bakta output: ${OUT}"
    rm -rf "${OUT}"

fi


# IMPORTANT:
# Do NOT create ${OUT}.
# Bakta creates its own output directory.

bakta \
        --db "${BAKTA_DB}" \
        --output "${OUT}" \
        --prefix "${ACC}" \
        --genus Rouxiella \
        --species badensis \
        --gram - \
        --threads "${THREADS}" \
        --keep-contig-headers \
        "${FNA}"


    # --------------------------------------------------------
    # Confirm required Bakta outputs exist
    # --------------------------------------------------------

    [[ -s "${OUT}/${ACC}.gff3" ]] || \
        die "Bakta GFF3 missing for ${ACC}"


    [[ -s "${OUT}/${ACC}.gbff" ]] || \
        die "Bakta GenBank file missing for ${ACC}"


    [[ -s "${OUT}/${ACC}.tsv" ]] || \
        die "Bakta TSV missing for ${ACC}"


    touch "${GENOME_DONE}"


    echo
    echo "${ACC}: annotation completed successfully."


done



# ------------------------------------------------------------
# Determine whether ALL genomes are finished
# ------------------------------------------------------------

ALL_ANNOTATED=true


for FNA in "${GENOME_FILES[@]}"; do

    ACC="$(basename "${FNA}" .fna)"

    if [[ ! -f "${BAKTA_ANNOT}/${ACC}/.annotation.done" ]]; then

        ALL_ANNOTATED=false

    fi

done


if [[ "${ALL_ANNOTATED}" == true ]]; then

    mark_done "05_bakta_annotation"

fi



# ============================================================
# STEP 06
# ANNOTATION SUMMARY
# ============================================================
#
# Generates:
#
#   all_feature_type_counts.tsv
#   annotation_feature_summary.tsv
#
# and separate extracted GFF files containing:
#
#   CDS
#   tRNA
#   ncRNA/sRNA
#   rRNA
#
# ============================================================

if ! step_done "06_annotation_summary"; then

    echo
    echo "============================================================"
    echo "STEP 06: Summarizing annotation features"
    echo "============================================================"
    echo


    GFF_FILES=( "${BAKTA_ANNOT}"/*/*.gff3 )


    if [[ "${#GFF_FILES[@]}" -eq 0 ]]; then

        die "No Bakta GFF3 files found."

    fi


    # --------------------------------------------------------
    # Detailed feature-type count table
    # --------------------------------------------------------

    ALL_TYPES="${FEATURE_RESULTS}/all_feature_type_counts.tsv"


    echo -e \
        "assembly\tfeature_type\tcount" \
        > "${ALL_TYPES}"


    # --------------------------------------------------------
    # Main requested-feature summary
    # --------------------------------------------------------

    SUMMARY="${FEATURE_RESULTS}/annotation_feature_summary.tsv"


    echo -e \
"assembly\tCDS\ttRNA\ttmRNA\trRNA\tncRNA\tncRNA_region\tCRISPR" \
        > "${SUMMARY}"


    for GFF in "${GFF_FILES[@]}"; do


        ACC="$(basename "${GFF}" .gff3)"


        # ----------------------------------------------------
        # Count every feature type Bakta produced
        # ----------------------------------------------------

        awk -F '\t' '

            BEGIN {
                OFS="\t"
            }

            !/^#/ && NF >= 3 {
                count[$3]++
            }

            END {

                for (type in count) {
                    print type, count[type]
                }

            }

        ' "${GFF}" \
            | sort -k1,1 \
            | while IFS=$'\t' read -r TYPE COUNT; do

                echo -e \
                    "${ACC}\t${TYPE}\t${COUNT}" \
                    >> "${ALL_TYPES}"

            done


        # ----------------------------------------------------
        # Feature counts of particular interest
        # ----------------------------------------------------

        CDS=$(awk -F '\t' \
            '$3=="CDS"{n++} END{print n+0}' \
            "${GFF}")


        TRNA=$(awk -F '\t' \
            '$3=="tRNA"{n++} END{print n+0}' \
            "${GFF}")


        TMRNA=$(awk -F '\t' \
            '$3=="tmRNA"{n++} END{print n+0}' \
            "${GFF}")


        RRNA=$(awk -F '\t' \
            '$3=="rRNA"{n++} END{print n+0}' \
            "${GFF}")


        NCRNA=$(awk -F '\t' \
            '$3=="ncRNA"{n++} END{print n+0}' \
            "${GFF}")


        NCRREG=$(awk -F '\t' '

            $3=="ncRNA-region" ||
            $3=="ncRNA_region" {
                n++
            }

            END {
                print n+0
            }

        ' "${GFF}")


        CRISPR=$(awk -F '\t' \
            '$3=="CRISPR"{n++} END{print n+0}' \
            "${GFF}")


        echo -e \
"${ACC}\t${CDS}\t${TRNA}\t${TMRNA}\t${RRNA}\t${NCRNA}\t${NCRREG}\t${CRISPR}" \
            >> "${SUMMARY}"


        # ----------------------------------------------------
        # Extract CDS features
        # ----------------------------------------------------

        awk -F '\t' '

            BEGIN {
                OFS="\t"
            }

            !/^#/ && $3=="CDS"

        ' "${GFF}" \
            > "${FEATURE_RESULTS}/${ACC}_CDS.gff3"


        # ----------------------------------------------------
        # Extract tRNA features
        # ----------------------------------------------------

        awk -F '\t' '

            BEGIN {
                OFS="\t"
            }

            !/^#/ && $3=="tRNA"

        ' "${GFF}" \
            > "${FEATURE_RESULTS}/${ACC}_tRNA.gff3"


        # ----------------------------------------------------
        # Extract rRNA features
        # ----------------------------------------------------

        awk -F '\t' '

            BEGIN {
                OFS="\t"
            }

            !/^#/ && $3=="rRNA"

        ' "${GFF}" \
            > "${FEATURE_RESULTS}/${ACC}_rRNA.gff3"


        # ----------------------------------------------------
        # Extract ncRNA / sRNA / regulatory RNA features
        # ----------------------------------------------------

        awk -F '\t' '

            BEGIN {
                OFS="\t"
            }

            !/^#/ &&
            (
                $3=="ncRNA" ||
                $3=="ncRNA-region" ||
                $3=="ncRNA_region" ||
                $3=="tmRNA"
            )

        ' "${GFF}" \
            > "${FEATURE_RESULTS}/${ACC}_ncRNA.gff3"


    done



    # ========================================================
    # INITIAL REGULATORY KEYWORD SCREEN
    # ========================================================

    CAND="${RESULTS}/regulatory_keyword_hits.txt"


    : > "${CAND}"


    TSV_FILES=( "${BAKTA_ANNOT}"/*/*.tsv )


    for TSV in "${TSV_FILES[@]}"; do


        ACC="$(basename "${TSV}" .tsv)"


        echo \
"==================== ${ACC} ====================" \
            >> "${CAND}"


        grep -iE \
'ArcZ|SlyA|LuxI|LuxR|ExpI|ExpR|Vfm|homoserine lactone|quorum' \
            "${TSV}" \
            >> "${CAND}" \
            || true


        echo >> "${CAND}"


    done


    echo
    echo "Annotation summary:"
    echo

    column -t "${SUMMARY}" || cat "${SUMMARY}"


    mark_done "06_annotation_summary"

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi



# ============================================================
# STEP 07
# ArcZ SEARCH USING RFAM / INFERNAL
# ============================================================
#
# ArcZ:
#
#   Rfam family RF00081
#
# Two searches are performed:
#
#   A. TRUSTED search
#      Uses the Rfam family gathering threshold (--cut_ga)
#
#   B. EXPLORATORY search
#      Reports weaker candidates down to E <= 1e-3
#
# This is much more appropriate for ArcZ than simply searching
# for a short nucleotide match with BLAST.
#
# ============================================================

if ! step_done "07_ArcZ_rfam"; then

    echo
    echo "============================================================"
    echo "STEP 07: ArcZ Rfam / Infernal search"
    echo "============================================================"
    echo


    ARCZ_CM="${ARCS}/ArcZ_RF00081.cm"


    # --------------------------------------------------------
    # Download ArcZ covariance model
    # --------------------------------------------------------

    if [[ ! -s "${ARCZ_CM}" ]]; then

        echo "Downloading Rfam ArcZ model RF00081"

        wget \
            -O "${ARCZ_CM}" \
            "https://rfam.org/family/RF00081/cm"

    fi


    [[ -s "${ARCZ_CM}" ]] || \
        die "ArcZ RF00081 covariance model was not downloaded."


    # --------------------------------------------------------
    # Search each genome independently
    # --------------------------------------------------------

    for FNA in "${GENOME_FILES[@]}"; do


        ACC="$(basename "${FNA}" .fna)"


        ARCZ_DONE="${ARCS}/${ACC}.done"


        if [[ -f "${ARCZ_DONE}" ]]; then

            echo
            echo "${ACC}: ArcZ search already completed — skipping"

            continue

        fi


        echo
        echo "Searching ArcZ in ${ACC}"


        # ====================================================
        # Trusted Rfam threshold
        # ====================================================

        cmsearch \
            --cpu "${THREADS}" \
            --cut_ga \
            --tblout "${ARCS}/${ACC}_ArcZ_trusted.tbl" \
            "${ARCZ_CM}" \
            "${FNA}" \
            > "${ARCS}/${ACC}_ArcZ_trusted.full.txt"


        # ====================================================
        # Exploratory search for weaker/distant candidates
        # ====================================================

        cmsearch \
            --cpu "${THREADS}" \
            -E 1e-3 \
            --tblout "${ARCS}/${ACC}_ArcZ_exploratory.tbl" \
            "${ARCZ_CM}" \
            "${FNA}" \
            > "${ARCS}/${ACC}_ArcZ_exploratory.full.txt"


        touch "${ARCZ_DONE}"


    done



    # ========================================================
    # BUILD COMBINED ArcZ SUMMARY
    # ========================================================

    HIT_COUNTS="${ARCS}/ArcZ_hit_counts.tsv"


    echo -e \
        "assembly\ttrusted_hits\texploratory_hits" \
        > "${HIT_COUNTS}"


    TRUSTED_SUMMARY="${ARCS}/ArcZ_trusted_all_genomes.txt"

    EXPLORATORY_SUMMARY="${ARCS}/ArcZ_exploratory_all_genomes.txt"


    : > "${TRUSTED_SUMMARY}"

    : > "${EXPLORATORY_SUMMARY}"


    for FNA in "${GENOME_FILES[@]}"; do


        ACC="$(basename "${FNA}" .fna)"


        TRUSTED_TBL="${ARCS}/${ACC}_ArcZ_trusted.tbl"

        EXPLORATORY_TBL="${ARCS}/${ACC}_ArcZ_exploratory.tbl"


        TRUSTED_COUNT=$(

            awk '

                !/^#/ && NF > 0 {
                    n++
                }

                END {
                    print n+0
                }

            ' "${TRUSTED_TBL}"

        )


        EXPLORATORY_COUNT=$(

            awk '

                !/^#/ && NF > 0 {
                    n++
                }

                END {
                    print n+0
                }

            ' "${EXPLORATORY_TBL}"

        )


        echo -e \
"${ACC}\t${TRUSTED_COUNT}\t${EXPLORATORY_COUNT}" \
            >> "${HIT_COUNTS}"


        echo \
"==================== ${ACC} ====================" \
            >> "${TRUSTED_SUMMARY}"


        grep -v '^#' \
            "${TRUSTED_TBL}" \
            >> "${TRUSTED_SUMMARY}" \
            || true


        echo \
            >> "${TRUSTED_SUMMARY}"


        echo \
"==================== ${ACC} ====================" \
            >> "${EXPLORATORY_SUMMARY}"


        grep -v '^#' \
            "${EXPLORATORY_TBL}" \
            >> "${EXPLORATORY_SUMMARY}" \
            || true


        echo \
            >> "${EXPLORATORY_SUMMARY}"


    done


    echo
    echo "ArcZ hit counts:"
    echo

    column -t "${HIT_COUNTS}" || cat "${HIT_COUNTS}"


    mark_done "07_ArcZ_rfam"

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
echo " ROUXIELLA PIPELINE FINISHED"
echo "============================================================"
echo

echo "Finished:"
echo "  $(date)"
echo

echo "Project:"
echo "  ${ROOT}"
echo

echo "NCBI genome table:"
echo "  ${NCBI}/R_badensis_assemblies.tsv"
echo

echo "20GA0316 NCBI match:"
echo "  ${NCBI}/R_badensis_20GA0316.tsv"
echo

echo "Accession list:"
echo "  ${NCBI}/accessions.txt"
echo

echo "Genome FASTA:"
echo "  ${GENOMES}/"
echo

echo "Original NCBI annotations:"
echo "  ${NCBI_ORIGINAL}/"
echo

echo "Bakta annotations:"
echo "  ${BAKTA_ANNOT}/"
echo

echo "Genome statistics:"
echo "  ${GENOME_QC}/genome_stats.tsv"
echo

echo "Annotation summary:"
echo "  ${FEATURE_RESULTS}/annotation_feature_summary.tsv"
echo

echo "All annotation feature types:"
echo "  ${FEATURE_RESULTS}/all_feature_type_counts.tsv"
echo

echo "Regulatory keyword hits:"
echo "  ${RESULTS}/regulatory_keyword_hits.txt"
echo

echo "ArcZ hit counts:"
echo "  ${ARCS}/ArcZ_hit_counts.tsv"
echo

echo "Trusted ArcZ hits:"
echo "  ${ARCS}/ArcZ_trusted_all_genomes.txt"
echo

echo "Exploratory ArcZ hits:"
echo "  ${ARCS}/ArcZ_exploratory_all_genomes.txt"
echo


echo "============================================================"
echo " FUTURE STEPS CAN BE ADDED BELOW THIS POINT"
echo "============================================================"



# ============================================================
# TEMPLATE FOR FUTURE ANALYSES
# ============================================================
#
# Example:
#
#
# if ! step_done "08_SlyA"; then
#
#     echo
#     echo "STEP 08: SlyA analysis"
#
#
#     YOUR_COMMANDS_HERE
#
#
#     mark_done "08_SlyA"
#
# else
#
#     echo
#     echo "STEP 08 already completed — skipping"
#
# fi
#
#
# ============================================================