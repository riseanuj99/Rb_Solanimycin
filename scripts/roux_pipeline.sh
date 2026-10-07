#!/usr/bin/env bash

#SBATCH --job-name=Roux_annotation
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail

# ============================================================
# ROUXIELLA BADENSIS GENOME / ANNOTATION PIPELINE
#
# Project:
#   /scratch/al98750/Roux
#
# Design:
#   Every completed stage gets a .done marker.
#   Re-running this script will SKIP completed stages.
#
# Current stages:
#   00 - software environment
#   01 - query current NCBI R. badensis assemblies
#   02 - download genomes + NCBI annotations
#   03 - prepare genome FASTA files
#   04 - download Bakta database
#   05 - comprehensive Bakta annotation
#   06 - summarize annotations / search regulatory genes
#   07 - ArcZ-specific Rfam covariance-model search
#
# Future analyses can be appended at the bottom.
# ============================================================


# ============================================================
# USER SETTINGS
# ============================================================

ROOT="/scratch/al98750/Roux"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

ENV_DIR="${ROOT}/envs/roux"

STATE="${ROOT}/.state"

NCBI="${ROOT}/01_NCBI"
GENOMES="${ROOT}/02_genomes"
BAKTA_DB_ROOT="${ROOT}/db/bakta"
BAKTA_DB="${BAKTA_DB_ROOT}/db-light"

ANNOT="${ROOT}/03_annotation/bakta"

RESULTS="${ROOT}/04_results"

ARCS="${RESULTS}/ArcZ"

LOGS="${ROOT}/logs"


mkdir -p \
    "${STATE}" \
    "${NCBI}" \
    "${GENOMES}" \
    "${BAKTA_DB_ROOT}" \
    "${ANNOT}" \
    "${RESULTS}" \
    "${ARCS}" \
    "${LOGS}"


echo
echo "============================================================"
echo " ROUXIELLA BADENSIS PIPELINE"
echo " Project: ${ROOT}"
echo " Threads: ${THREADS}"
echo " Started: $(date)"
echo "============================================================"
echo


# ============================================================
# HELPER FUNCTIONS
# ============================================================

done_step () {

    [[ -f "${STATE}/$1.done" ]]

}


mark_done () {

    date > "${STATE}/$1.done"

}


# ============================================================
# STEP 00
# CREATE SOFTWARE ENVIRONMENT
#
# Software:
#   NCBI datasets
#   Bakta
#   BLAST+
#   Infernal
#   SeqKit
#   jq
# ============================================================

if ! done_step "00_environment"; then

    echo ">>> STEP 00: Creating software environment"

    source "$(conda info --base)/etc/profile.d/conda.sh"

    if [[ ! -d "${ENV_DIR}" ]]; then

        if command -v mamba >/dev/null 2>&1; then
            SOLVER="mamba"
        else
            SOLVER="conda"
        fi

        "${SOLVER}" create -y \
            -p "${ENV_DIR}" \
            -c conda-forge \
            -c bioconda \
            ncbi-datasets-cli \
            bakta \
            blast \
            infernal \
            seqkit \
            jq \
            wget

    fi

    mark_done "00_environment"

else

    echo ">>> STEP 00 already completed — skipping"

fi


source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_DIR}"


echo
echo "Software versions:"
datasets version || true
bakta --version || true
blastn -version | head -1 || true
cmsearch -h | head -2 || true
echo


# ============================================================
# STEP 01
# QUERY NCBI FOR ALL CURRENT ROUXIELLA BADENSIS ASSEMBLIES
#
# IMPORTANT:
# We do NOT hard-code only the genomes known today.
# NCBI is queried for the exact taxon "Rouxiella badensis".
# ============================================================

if ! done_step "01_ncbi_query"; then

    echo ">>> STEP 01: Querying NCBI for Rouxiella badensis genomes"

    datasets summary genome taxon \
        "Rouxiella badensis" \
        --tax-exact-match \
        --as-json-lines \
        > "${NCBI}/R_badensis_metadata.jsonl"


    datasets summary genome taxon \
        "Rouxiella badensis" \
        --tax-exact-match \
        --as-json-lines \
    | dataformat tsv genome \
        --fields accession,organism-name,assminfo-name \
        > "${NCBI}/R_badensis_assemblies.tsv"


    # Extract assembly accessions.
    tail -n +2 "${NCBI}/R_badensis_assemblies.tsv" \
        | cut -f1 \
        | grep -E '^GC[AF]_' \
        | sort -u \
        > "${NCBI}/accessions.txt"


    echo
    echo "Assemblies found:"
    cat "${NCBI}/R_badensis_assemblies.tsv"
    echo

    echo "Total assemblies:"
    wc -l "${NCBI}/accessions.txt"

    mark_done "01_ncbi_query"

else

    echo ">>> STEP 01 already completed — skipping"

fi


# ============================================================
# STEP 02
# DOWNLOAD GENOMES + ORIGINAL NCBI ANNOTATIONS
#
# We retain:
#   genome FASTA
#   GFF3
#   GenBank file
#   RNA FASTA
#   CDS FASTA
#   protein FASTA
#
# Thus we preserve NCBI/PGAP annotation in addition to our
# standardized Bakta reannotation.
# ============================================================

if ! done_step "02_ncbi_download"; then

    echo ">>> STEP 02: Downloading NCBI genomes and annotations"

    rm -f "${NCBI}/R_badensis_ncbi.zip"

    datasets download genome accession \
        --inputfile "${NCBI}/accessions.txt" \
        --include genome,gff3,gbff,rna,cds,protein \
        --filename "${NCBI}/R_badensis_ncbi.zip"


    rm -rf "${NCBI}/package"

    mkdir -p "${NCBI}/package"

    unzip -q \
        "${NCBI}/R_badensis_ncbi.zip" \
        -d "${NCBI}/package"


    mark_done "02_ncbi_download"

else

    echo ">>> STEP 02 already completed — skipping"

fi


# ============================================================
# STEP 03
# PREPARE ONE CLEAN GENOME FASTA PER ASSEMBLY
# ============================================================

if ! done_step "03_prepare_genomes"; then

    echo ">>> STEP 03: Preparing genome FASTA files"

    while read -r ACC; do

        SOURCE_DIR="${NCBI}/package/ncbi_dataset/data/${ACC}"

        if [[ ! -d "${SOURCE_DIR}" ]]; then

            echo "WARNING: No NCBI directory for ${ACC}"
            continue

        fi


        FNA=$(find "${SOURCE_DIR}" \
            -maxdepth 1 \
            -type f \
            -name "*_genomic.fna" \
            | head -1)


        if [[ -z "${FNA}" ]]; then

            echo "WARNING: Genome FASTA not found for ${ACC}"
            continue

        fi


        ln -sfn \
            "${FNA}" \
            "${GENOMES}/${ACC}.fna"


        echo "${ACC} -> ${FNA}"

    done < "${NCBI}/accessions.txt"


    echo
    echo "Prepared genomes:"
    ls -lh "${GENOMES}"
    echo


    # Basic genome statistics
    seqkit stats \
        "${GENOMES}"/*.fna \
        > "${RESULTS}/genome_stats.tsv"


    mark_done "03_prepare_genomes"

else

    echo ">>> STEP 03 already completed — skipping"

fi


# ============================================================
# STEP 04
# DOWNLOAD BAKTA DATABASE
#
# Light database:
#   ~4 GB unpacked
#
# IMPORTANT:
# Light DB still contains ALL non-coding RNA information.
# It is sufficient for:
#   sRNA
#   ncRNA
#   tRNA
#   tmRNA
#   rRNA
#   CDS
#   pseudogene
#   sORF
#   CRISPR
#
# Full Bakta DB can be used later if we want more detailed
# protein functional annotation.
# ============================================================

if [[ ! -d "${BAKTA_DB}" ]]; then

    echo ">>> STEP 04: Downloading Bakta LIGHT database"

    bakta_db download \
        --output "${BAKTA_DB_ROOT}" \
        --type light

else

    echo ">>> STEP 04 Bakta database already present — skipping"

fi


# ============================================================
# STEP 05
# COMPREHENSIVE ANNOTATION WITH BAKTA
#
# Includes:
#   CDS
#   small ORFs
#   pseudogenes
#   tRNA
#   tmRNA
#   rRNA
#   ncRNA / sRNA
#   antisense RNAs
#   ribozymes
#   antitoxins
#   riboswitches
#   RNA leaders
#   CRISPR
#   oriC / oriT
#
# Each genome gets its OWN completion marker.
# ============================================================

echo
echo ">>> STEP 05: Bakta annotations"
echo


for FNA in "${GENOMES}"/*.fna; do

    ACC=$(basename "${FNA}" .fna)

    OUT="${ANNOT}/${ACC}"

    GENOME_DONE="${OUT}/.done"


    if [[ -f "${GENOME_DONE}" ]]; then

        echo "${ACC}: annotation already completed — skipping"
        continue

    fi


    echo
    echo "------------------------------------------------------------"
    echo "Annotating ${ACC}"
    echo "------------------------------------------------------------"


    # Remove only an incomplete previous annotation.
    if [[ -d "${OUT}" ]]; then
        rm -rf "${OUT}"
    fi


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


    touch "${GENOME_DONE}"

done


mark_done "05_bakta_annotation"


# ============================================================
# STEP 06
# CREATE ANNOTATION SUMMARY AND LOOK FOR GENES OF INTEREST
# ============================================================

if ! done_step "06_annotation_summary"; then

    echo
    echo ">>> STEP 06: Summarizing annotation"
    echo


    SUMMARY="${RESULTS}/annotation_feature_counts.tsv"

    echo -e \
        "assembly\tCDS\tsORF\ttRNA\ttmRNA\trRNA\tncRNA\tncRNA_region\tCRISPR" \
        > "${SUMMARY}"


    for GFF in "${ANNOT}"/*/*.gff3; do

        ACC=$(basename "${GFF}" .gff3)


        CDS=$(awk -F'\t' '$3=="CDS"{n++} END{print n+0}' "${GFF}")
        SORF=$(awk -F'\t' '$3=="sORF"{n++} END{print n+0}' "${GFF}")
        TRNA=$(awk -F'\t' '$3=="tRNA"{n++} END{print n+0}' "${GFF}")
        TMRNA=$(awk -F'\t' '$3=="tmRNA"{n++} END{print n+0}' "${GFF}")
        RRNA=$(awk -F'\t' '$3=="rRNA"{n++} END{print n+0}' "${GFF}")
        NCRNA=$(awk -F'\t' '$3=="ncRNA"{n++} END{print n+0}' "${GFF}")
        NCRREG=$(awk -F'\t' '$3=="ncRNA-region"{n++} END{print n+0}' "${GFF}")
        CRISPR=$(awk -F'\t' '$3=="CRISPR"{n++} END{print n+0}' "${GFF}")


        echo -e \
            "${ACC}\t${CDS}\t${SORF}\t${TRNA}\t${TMRNA}\t${RRNA}\t${NCRNA}\t${NCRREG}\t${CRISPR}" \
            >> "${SUMMARY}"

    done


    # --------------------------------------------------------
    # Initial keyword search.
    #
    # This is only a FIRST SCREEN.
    # Later we will do sequence/homology-based searches for
    # LuxI/R, Vfm, SlyA, sol, etc.
    # --------------------------------------------------------

    CAND="${RESULTS}/regulatory_keyword_hits.txt"

    : > "${CAND}"


    for TSV in "${ANNOT}"/*/*.tsv; do

        ACC=$(basename "${TSV}" .tsv)

        echo \
            "==================== ${ACC} ====================" \
            >> "${CAND}"


        grep -iE \
            'ArcZ|SlyA|LuxI|LuxR|ExpI|ExpR|Vfm|homoserine lactone' \
            "${TSV}" \
            >> "${CAND}" \
            || true


        echo >> "${CAND}"

    done


    mark_done "06_annotation_summary"

else

    echo ">>> STEP 06 already completed — skipping"

fi


# ============================================================
# STEP 07
# ARcZ-SPECIFIC RFAM SEARCH
#
# ArcZ = Rfam family RF00081
#
# This is much more appropriate than a short BLAST hit because
# Infernal uses the conserved RNA sequence + secondary-structure
# covariance model.
# ============================================================

if ! done_step "07_ArcZ_rfam"; then

    echo
    echo ">>> STEP 07: ArcZ-specific Rfam/Infernal search"
    echo


    ARCZ_CM="${ARCS}/ArcZ_RF00081.cm"


    if [[ ! -s "${ARCZ_CM}" ]]; then

        wget \
            -O "${ARCZ_CM}" \
            "https://rfam.org/family/RF00081/cm"

    fi


    COMBINED="${ARCS}/ArcZ_all_genomes_summary.txt"

    : > "${COMBINED}"


    for FNA in "${GENOMES}"/*.fna; do

        ACC=$(basename "${FNA}" .fna)


        echo
        echo "Searching ArcZ in ${ACC}"


        cmsearch \
            --cut_ga \
            --tblout "${ARCS}/${ACC}_ArcZ.tbl" \
            "${ARCZ_CM}" \
            "${FNA}" \
            > "${ARCS}/${ACC}_ArcZ.full.txt"


        echo \
            "==================== ${ACC} ====================" \
            >> "${COMBINED}"


        grep -v '^#' \
            "${ARCS}/${ACC}_ArcZ.tbl" \
            >> "${COMBINED}" \
            || true


        echo >> "${COMBINED}"

    done


    mark_done "07_ArcZ_rfam"

else

    echo ">>> STEP 07 already completed — skipping"

fi


# ============================================================
# FINAL REPORT
# ============================================================

echo
echo "============================================================"
echo " PIPELINE COMPLETE"
echo "============================================================"
echo
echo "NCBI assembly table:"
echo "  ${NCBI}/R_badensis_assemblies.tsv"
echo
echo "Genome FASTAs:"
echo "  ${GENOMES}/"
echo
echo "Bakta annotation:"
echo "  ${ANNOT}/"
echo
echo "Feature count table:"
echo "  ${RESULTS}/annotation_feature_counts.tsv"
echo
echo "Initial regulatory keyword search:"
echo "  ${RESULTS}/regulatory_keyword_hits.txt"
echo
echo "ArcZ covariance-model results:"
echo "  ${ARCS}/ArcZ_all_genomes_summary.txt"
echo
echo "Finished: $(date)"
echo


# ============================================================
# ADD FUTURE STEPS BELOW THIS LINE
#
# Template:
#
# if ! done_step "08_example"; then
#
#     echo ">>> STEP 08"
#
#     YOUR_COMMANDS_HERE
#
#     mark_done "08_example"
#
# else
#
#     echo ">>> STEP 08 already completed — skipping"
#
# fi
#
# ============================================================