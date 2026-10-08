#!/usr/bin/env bash

#SBATCH --job-name=Roux_Vfm
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/Roux_Vfm_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/Roux_Vfm_%j.err

set -euo pipefail
shopt -s nullglob


# ============================================================
# OBJECTIVE 2
# VFM QUORUM-SENSING SYSTEM IN ROUXIELLA BADENSIS
# ============================================================
#
# Question:
#
# Has R. badensis retained components of the Vfm regulatory
# system that regulates solanimycin in Dickeya?
#
#
# Strategy:
#
#   Canonical Vfm reference:
#       Dickeya dadantii 3937
#       NC_014500.1
#
#   Positive control:
#       Dickeya solani MK10
#       GCA_000365285.1
#
#   Primary target:
#       Rouxiella badensis 20GA0316
#       GCF_020740305.1
#
#   Also examine all unique R. badensis genomes currently
#   present in:
#
#       /scratch/al98750/Roux/02_genomes/
#
#
# PRIMARY REGULATORY COMPONENTS:
#
#   VfmI = sensor histidine kinase
#   VfmH = response regulator
#   VfmE = AraC-family transcriptional regulator
#
#
# We test:
#
#   1. sequence homology
#   2. query coverage
#   3. genomic coordinates
#   4. whether VfmE/H/I remain clustered
#   5. conservation of the surrounding canonical Vfm region
#
#
# CHECKPOINT SYSTEM:
#
# Completed steps are marked under:
#
#   /scratch/al98750/Roux/06_QS_Vfm/.state/
#
# Rerunning this script skips finished steps.
#
# Future steps can be appended as STEP 08, STEP 09, etc.
#
# ============================================================



# ============================================================
# PATHS
# ============================================================

ROOT="/scratch/al98750/Roux"

WORK="${ROOT}/06_QS_Vfm"

REFERENCE="${WORK}/01_reference"

DDAD="${REFERENCE}/Ddadantii_3937"

MK10="${REFERENCE}/Dsolani_MK10"

QUERY="${WORK}/02_queries"

DB="${WORK}/03_BLAST_databases"

RAW="${WORK}/04_raw_BLAST"

RESULTS="${WORK}/05_results"

CONTEXT="${WORK}/06_genomic_context"

STATE="${WORK}/.state"

TMP="${WORK}/tmp"


THREADS="${SLURM_CPUS_PER_TASK:-4}"


# ============================================================
# ACCESSIONS
# ============================================================

DDAD_NUCCORE="NC_014500.1"

MK10_ACC="GCA_000365285.1"

RB_PRIMARY="GCF_020740305.1"


# ============================================================
# CANONICAL D. DADANTII LOCUS TAGS
#
# Literature / current RefSeq annotation:
#
#   DDA3937_RS20815 = VfmI
#   DDA3937_RS20820 = VfmH
#   DDA3937_RS20835 = VfmE
#
# ============================================================

VFMI_TAG="DDA3937_RS20815"
VFMH_TAG="DDA3937_RS20820"
VFME_TAG="DDA3937_RS20835"


# ============================================================
# VFM REGION WINDOW
#
# The entire Vfm locus is about 25-30 kb.
#
# We identify the canonical VfmE/H/I coordinates and retain
# surrounding CDSs within +/- 20 kb.
#
# This gives us a Vfm-region homology/synteny screen in
# addition to the three regulatory proteins.
#
# ============================================================

VFM_REGION_PADDING=20000



# ============================================================
# CREATE ALL OUTPUT DIRECTORIES
# ============================================================

mkdir -p \
    "${WORK}" \
    "${REFERENCE}" \
    "${DDAD}" \
    "${MK10}" \
    "${QUERY}" \
    "${DB}" \
    "${RAW}" \
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
# ACTIVATE EXISTING ROUX ENVIRONMENT
# ============================================================

source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate "${ROOT}/envs/roux"



# ============================================================
# JOB INFORMATION
# ============================================================

echo
echo "============================================================"
echo " ROUXIELLA VFM ANALYSIS"
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
echo "Canonical Vfm reference:"
echo "  D. dadantii 3937 (${DDAD_NUCCORE})"

echo
echo "Positive control:"
echo "  D. solani MK10 (${MK10_ACC})"

echo
echo "Primary Rouxiella:"
echo "  ${RB_PRIMARY}"

echo
echo "Started:"
echo "  $(date)"

echo



# ============================================================
# STEP 01
# SOFTWARE + INPUT CHECK
# ============================================================

if ! step_done "01_environment"; then

    echo
    echo "============================================================"
    echo "STEP 01: Environment checks"
    echo "============================================================"
    echo


    for PROGRAM in \
        datasets \
        makeblastdb \
        tblastn \
        wget \
        unzip \
        python

    do

        command -v "${PROGRAM}" >/dev/null 2>&1 || \
            die "${PROGRAM} was not found."

        echo "FOUND: ${PROGRAM}"

    done


    RB_GENOMES=( "${ROOT}/02_genomes/"*.fna )


    if [[ "${#RB_GENOMES[@]}" -eq 0 ]]; then

        die "No R. badensis genomes found in ${ROOT}/02_genomes"

    fi


    [[ -s "${ROOT}/02_genomes/${RB_PRIMARY}.fna" ]] || \
        die "Primary 20GA0316 genome ${RB_PRIMARY} is missing."


    echo
    echo "R. badensis genomes available:"
    echo "  ${#RB_GENOMES[@]}"
    echo


    mark_done "01_environment"

else

    echo
    echo ">>> STEP 01 already completed — skipping"

fi



# ============================================================
# STEP 02
# DOWNLOAD CANONICAL D. DADANTII 3937 CDS PROTEINS
# ============================================================

if ! step_done "02_Ddadantii_reference"; then

    echo
    echo "============================================================"
    echo "STEP 02: Downloading canonical D. dadantii Vfm reference"
    echo "============================================================"
    echo


    DDAD_CDS="${DDAD}/NC_014500.1_cds_proteins.faa"


    wget -q \
        -O "${DDAD_CDS}" \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${DDAD_NUCCORE}&rettype=fasta_cds_aa&retmode=text"


    [[ -s "${DDAD_CDS}" ]] || \
        die "D. dadantii CDS protein download failed."


    grep -q '^>' "${DDAD_CDS}" || \
        die "D. dadantii protein FASTA is invalid."


    echo "Downloaded:"
    echo "  ${DDAD_CDS}"

    echo
    echo "Number of proteins:"

    grep -c '^>' "${DDAD_CDS}"


    mark_done "02_Ddadantii_reference"

else

    echo
    echo ">>> STEP 02 already completed — skipping"

fi



DDAD_CDS="${DDAD}/NC_014500.1_cds_proteins.faa"



# ============================================================
# STEP 03
# EXTRACT VfmE / VfmH / VfmI AND CANONICAL VFM REGION
# ============================================================

if ! step_done "03_extract_Vfm_reference"; then

    echo
    echo "============================================================"
    echo "STEP 03: Extracting canonical Vfm reference proteins"
    echo "============================================================"
    echo


    export DDAD_CDS
    export QUERY
    export RESULTS
    export VFM_REGION_PADDING

    python <<'PY'

from pathlib import Path
import os
import re
import sys


fasta = Path(os.environ["DDAD_CDS"])

query_dir = Path(os.environ["QUERY"])

results_dir = Path(os.environ["RESULTS"])

padding = int(os.environ["VFM_REGION_PADDING"])


# ------------------------------------------------------------
# Known regulatory anchors
# ------------------------------------------------------------

anchor_tags = {

    "vfmI": "DDA3937_RS20815",

    "vfmH": "DDA3937_RS20820",

    "vfmE": "DDA3937_RS20835",

}


# ------------------------------------------------------------
# FASTA parser
# ------------------------------------------------------------

def read_fasta(path):

    records = []

    header = None

    seq = []


    with open(path) as handle:

        for line in handle:

            line = line.rstrip()


            if line.startswith(">"):

                if header is not None:

                    records.append(
                        (header, "".join(seq))
                    )


                header = line[1:]

                seq = []


            else:

                seq.append(line.strip())


        if header is not None:

            records.append(
                (header, "".join(seq))
            )


    return records



records = read_fasta(fasta)


def get_tag(header):

    m = re.search(
        r"\[locus_tag=([^\]]+)\]",
        header
    )

    return m.group(1) if m else None



def get_coords(header):

    m = re.search(
        r"\[location=([^\]]+)\]",
        header
    )


    if not m:

        return None


    location = m.group(1)


    nums = [
        int(x)
        for x in re.findall(r"\d+", location)
    ]


    if len(nums) < 2:

        return None


    return min(nums), max(nums)



# ------------------------------------------------------------
# Find VfmE/H/I
# ------------------------------------------------------------

anchors = {}


for header, seq in records:

    tag = get_tag(header)


    for gene, expected_tag in anchor_tags.items():

        if (
            tag == expected_tag
            or expected_tag.lower() in header.lower()
            or f"gene={gene}".lower() in header.lower()
            or gene.lower() in header.lower()
        ):

            coords = get_coords(header)


            anchors[gene] = {

                "tag": tag or expected_tag,

                "header": header,

                "sequence": seq,

                "coords": coords,

            }



missing = [
    gene
    for gene in anchor_tags
    if gene not in anchors
]


if missing:

    print(
        "ERROR: missing canonical Vfm anchors:",
        ", ".join(missing),
        file=sys.stderr,
    )

    sys.exit(1)



# ------------------------------------------------------------
# Write clean VfmE/H/I query FASTA
# ------------------------------------------------------------

key_fasta = query_dir / "Ddadantii_VfmEHI.faa"


with open(key_fasta, "w") as out:

    for gene in ("vfmE", "vfmH", "vfmI"):

        rec = anchors[gene]

        out.write(
            f">{gene}|{rec['tag']}|Ddadantii3937\n"
        )


        sequence = rec["sequence"]


        for i in range(0, len(sequence), 70):

            out.write(
                sequence[i:i+70] + "\n"
            )



# ------------------------------------------------------------
# Anchor information
# ------------------------------------------------------------

anchor_tsv = results_dir / "Ddadantii_Vfm_key_reference.tsv"


with open(anchor_tsv, "w") as out:

    out.write(
        "gene\tlocus_tag\tstart\tend\tprotein_length\n"
    )


    for gene in ("vfmE", "vfmH", "vfmI"):

        rec = anchors[gene]

        start, end = rec["coords"]


        out.write(
            f"{gene}\t"
            f"{rec['tag']}\t"
            f"{start}\t"
            f"{end}\t"
            f"{len(rec['sequence'])}\n"
        )



# ------------------------------------------------------------
# Define canonical Vfm-region window
# ------------------------------------------------------------

starts = [
    rec["coords"][0]
    for rec in anchors.values()
]

ends = [
    rec["coords"][1]
    for rec in anchors.values()
]


region_start = max(
    1,
    min(starts) - padding
)

region_end = max(ends) + padding


print(
    f"Canonical Vfm-region window: "
    f"{region_start}-{region_end}"
)



# ------------------------------------------------------------
# Extract all CDS proteins within this region
#
# This provides a broad Vfm-locus conservation/synteny screen.
# The E/H/I analysis remains the primary regulatory test.
# ------------------------------------------------------------

region_records = []


for header, seq in records:

    coords = get_coords(header)


    if coords is None:

        continue


    start, end = coords


    if (
        start <= region_end
        and end >= region_start
    ):

        region_records.append(
            (
                header,
                seq,
                start,
                end,
                get_tag(header),
            )
        )



region_records.sort(
    key=lambda x: x[2]
)


region_fasta = query_dir / "Ddadantii_Vfm_region_proteins.faa"

region_tsv = results_dir / "Ddadantii_Vfm_region_reference.tsv"


with open(region_fasta, "w") as fa, \
     open(region_tsv, "w") as table:


    table.write(
        "query_id\tlocus_tag\tstart\tend\tprotein_length\n"
    )


    for i, (
        header,
        seq,
        start,
        end,
        tag,
    ) in enumerate(region_records, start=1):


        query_id = (
            tag
            if tag
            else f"VFMregion_{i:02d}"
        )


        fa.write(
            f">{query_id}\n"
        )


        for j in range(0, len(seq), 70):

            fa.write(
                seq[j:j+70] + "\n"
            )


        table.write(
            f"{query_id}\t"
            f"{tag or 'NA'}\t"
            f"{start}\t"
            f"{end}\t"
            f"{len(seq)}\n"
        )



with open(
    results_dir / "Ddadantii_Vfm_region_coordinates.txt",
    "w",
) as out:

    out.write(
        f"{region_start}\t{region_end}\n"
    )


print(
    f"Extracted {len(region_records)} "
    f"proteins from canonical Vfm region."
)

PY


    [[ -s "${QUERY}/Ddadantii_VfmEHI.faa" ]] || \
        die "VfmE/H/I queries were not created."


    [[ -s "${QUERY}/Ddadantii_Vfm_region_proteins.faa" ]] || \
        die "Canonical Vfm-region protein set was not created."


    echo
    echo "Key Vfm reference:"
    cat "${RESULTS}/Ddadantii_Vfm_key_reference.tsv"

    echo


    mark_done "03_extract_Vfm_reference"

else

    echo
    echo ">>> STEP 03 already completed — skipping"

fi



# ============================================================
# STEP 04
# DOWNLOAD D. SOLANI MK10 POSITIVE-CONTROL GENOME
# ============================================================

if ! step_done "04_download_MK10"; then

    echo
    echo "============================================================"
    echo "STEP 04: Downloading D. solani MK10 positive control"
    echo "============================================================"
    echo


    ZIP="${MK10}/MK10.zip"

    PACKAGE="${MK10}/package"


    rm -f "${ZIP}"

    rm -rf "${PACKAGE}"


    datasets download genome accession \
        "${MK10_ACC}" \
        --include genome \
        --filename "${ZIP}"


    [[ -s "${ZIP}" ]] || \
        die "MK10 genome download failed."


    mkdir -p "${PACKAGE}"


    unzip -q \
        "${ZIP}" \
        -d "${PACKAGE}"


    MK10_FNA=$(

        find "${PACKAGE}" \
            -type f \
            -name "*.fna" \
            | head -1

    )


    [[ -n "${MK10_FNA}" && -s "${MK10_FNA}" ]] || \
        die "MK10 genome FASTA was not found."


    cp -f \
        "${MK10_FNA}" \
        "${MK10}/Dsolani_MK10.fna"


    mark_done "04_download_MK10"

else

    echo
    echo ">>> STEP 04 already completed — skipping"

fi



MK10_FNA="${MK10}/Dsolani_MK10.fna"



# ============================================================
# STEP 05
# BUILD NUCLEOTIDE BLAST DATABASES
# ============================================================

if ! step_done "05_BLAST_databases"; then

    echo
    echo "============================================================"
    echo "STEP 05: Building nucleotide BLAST databases"
    echo "============================================================"
    echo


    # --------------------------------------------------------
    # D. solani MK10 positive control
    # --------------------------------------------------------

    mkdir -p "${DB}/Dsolani_MK10"


    makeblastdb \
        -in "${MK10_FNA}" \
        -dbtype nucl \
        -parse_seqids \
        -out "${DB}/Dsolani_MK10/MK10"


    # --------------------------------------------------------
    # All R. badensis genomes
    # --------------------------------------------------------

    RB_GENOMES=( "${ROOT}/02_genomes/"*.fna )


    for FNA in "${RB_GENOMES[@]}"; do


        ACC="$(basename "${FNA}" .fna)"


        echo "Building DB: ${ACC}"


        mkdir -p "${DB}/${ACC}"


        makeblastdb \
            -in "${FNA}" \
            -dbtype nucl \
            -parse_seqids \
            -out "${DB}/${ACC}/${ACC}"


    done


    mark_done "05_BLAST_databases"

else

    echo
    echo ">>> STEP 05 already completed — skipping"

fi



# ============================================================
# STEP 06
# SEARCH VfmE / VfmH / VfmI
# ============================================================

if ! step_done "06_key_Vfm_search"; then

    echo
    echo "============================================================"
    echo "STEP 06: Searching VfmE/H/I"
    echo "============================================================"
    echo


    KEY_QUERY="${QUERY}/Ddadantii_VfmEHI.faa"


    # --------------------------------------------------------
    # Positive control
    # --------------------------------------------------------

    tblastn \
        -query "${KEY_QUERY}" \
        -db "${DB}/Dsolani_MK10/MK10" \
        -evalue 1e-5 \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -num_threads "${THREADS}" \
        -outfmt \
'6 qseqid sseqid pident length qlen qcovhsp sstart send evalue bitscore' \
        > "${RAW}/key_Vfm_vs_Dsolani_MK10.tsv"


    # --------------------------------------------------------
    # R. badensis genomes
    # --------------------------------------------------------

    RB_GENOMES=( "${ROOT}/02_genomes/"*.fna )


    for FNA in "${RB_GENOMES[@]}"; do


        ACC="$(basename "${FNA}" .fna)"


        SEARCH_DONE="${STATE}/06_key_${ACC}.done"


        if [[ -f "${SEARCH_DONE}" ]]; then

            echo "${ACC}: key Vfm search already done — skipping"

            continue

        fi


        echo "Searching VfmE/H/I in ${ACC}"


        tblastn \
            -query "${KEY_QUERY}" \
            -db "${DB}/${ACC}/${ACC}" \
            -evalue 1e-5 \
            -max_target_seqs 10 \
            -max_hsps 1 \
            -num_threads "${THREADS}" \
            -outfmt \
'6 qseqid sseqid pident length qlen qcovhsp sstart send evalue bitscore' \
            > "${RAW}/key_Vfm_vs_${ACC}.tsv"


        touch "${SEARCH_DONE}"


    done


    mark_done "06_key_Vfm_search"

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi



# ============================================================
# STEP 07
# SEARCH THE BROADER CANONICAL VFM REGION
# ============================================================

if ! step_done "07_Vfm_region_search"; then

    echo
    echo "============================================================"
    echo "STEP 07: Searching broader canonical Vfm region"
    echo "============================================================"
    echo


    REGION_QUERY="${QUERY}/Ddadantii_Vfm_region_proteins.faa"


    # --------------------------------------------------------
    # Positive-control D. solani
    # --------------------------------------------------------

    tblastn \
        -query "${REGION_QUERY}" \
        -db "${DB}/Dsolani_MK10/MK10" \
        -evalue 1e-5 \
        -max_target_seqs 5 \
        -max_hsps 1 \
        -num_threads "${THREADS}" \
        -outfmt \
'6 qseqid sseqid pident length qlen qcovhsp sstart send evalue bitscore' \
        > "${RAW}/Vfm_region_vs_Dsolani_MK10.tsv"


    # --------------------------------------------------------
    # R. badensis
    # --------------------------------------------------------

    RB_GENOMES=( "${ROOT}/02_genomes/"*.fna )


    for FNA in "${RB_GENOMES[@]}"; do


        ACC="$(basename "${FNA}" .fna)"


        SEARCH_DONE="${STATE}/07_region_${ACC}.done"


        if [[ -f "${SEARCH_DONE}" ]]; then

            echo "${ACC}: Vfm-region search already done — skipping"

            continue

        fi


        echo "Searching Vfm region in ${ACC}"


        tblastn \
            -query "${REGION_QUERY}" \
            -db "${DB}/${ACC}/${ACC}" \
            -evalue 1e-5 \
            -max_target_seqs 5 \
            -max_hsps 1 \
            -num_threads "${THREADS}" \
            -outfmt \
'6 qseqid sseqid pident length qlen qcovhsp sstart send evalue bitscore' \
            > "${RAW}/Vfm_region_vs_${ACC}.tsv"


        touch "${SEARCH_DONE}"


    done


    mark_done "07_Vfm_region_search"

else

    echo
    echo ">>> STEP 07 already completed — skipping"

fi



# ============================================================
# STEP 08
# SUMMARIZE HOMOLOGY + GENOMIC ORGANIZATION
# ============================================================

if ! step_done "08_Vfm_summary"; then

    echo
    echo "============================================================"
    echo "STEP 08: Summarizing Vfm conservation and organization"
    echo "============================================================"
    echo


    export ROOT
    export RAW
    export RESULTS
    export RB_PRIMARY


    python <<'PY'

from pathlib import Path
from collections import defaultdict
import os


root = Path(os.environ["ROOT"])

raw = Path(os.environ["RAW"])

results = Path(os.environ["RESULTS"])

primary = os.environ["RB_PRIMARY"]


# ============================================================
# Genome list
# ============================================================

genomes = ["Dsolani_MK10"]

genomes += sorted(
    p.stem
    for p in (root / "02_genomes").glob("*.fna")
)



# ============================================================
# Thresholds
#
# These are screening thresholds, NOT proof of orthology.
#
# High confidence:
#   >= 40% identity
#   >= 70% query coverage
#   E <= 1e-20
#
# Candidate:
#   >= 30% identity
#   >= 50% query coverage
#   E <= 1e-10
#
# Final inference also depends heavily on synteny.
# ============================================================

def classify(identity, qcov, evalue):

    if (
        identity >= 40
        and qcov >= 70
        and evalue <= 1e-20
    ):

        return "HIGH"


    if (
        identity >= 30
        and qcov >= 50
        and evalue <= 1e-10
    ):

        return "CANDIDATE"


    return "WEAK"



def read_best_hits(path):

    best = {}


    if not path.exists():

        return best


    with open(path) as handle:

        for line in handle:

            if not line.strip():

                continue


            parts = line.rstrip().split("\t")


            (
                query,
                target,
                pident,
                length,
                qlen,
                qcov,
                sstart,
                send,
                evalue,
                bitscore,
            ) = parts


            rec = {

                "query": query,

                "target": target,

                "identity": float(pident),

                "length": int(length),

                "qlen": int(qlen),

                "qcov": float(qcov),

                "sstart": int(sstart),

                "send": int(send),

                "evalue": float(evalue),

                "bitscore": float(bitscore),

            }


            if (
                query not in best
                or rec["bitscore"]
                > best[query]["bitscore"]
            ):

                best[query] = rec


    return best



# ============================================================
# KEY VfmE/H/I RESULTS
# ============================================================

key_rows = []

architecture_rows = []


for genome in genomes:


    if genome == "Dsolani_MK10":

        path = raw / "key_Vfm_vs_Dsolani_MK10.tsv"

    else:

        path = raw / f"key_Vfm_vs_{genome}.tsv"


    hits = read_best_hits(path)


    accepted = []


    for gene in ("vfmE", "vfmH", "vfmI"):


        possible = [
            rec
            for query, rec in hits.items()
            if query.startswith(gene + "|")
        ]


        if possible:

            rec = max(
                possible,
                key=lambda x: x["bitscore"],
            )


            status = classify(
                rec["identity"],
                rec["qcov"],
                rec["evalue"],
            )


            start = min(
                rec["sstart"],
                rec["send"],
            )

            end = max(
                rec["sstart"],
                rec["send"],
            )

            strand = (
                "+"
                if rec["sstart"] <= rec["send"]
                else "-"
            )


            key_rows.append({

                "genome": genome,

                "gene": gene,

                "target": rec["target"],

                "identity": rec["identity"],

                "qcov": rec["qcov"],

                "start": start,

                "end": end,

                "strand": strand,

                "evalue": rec["evalue"],

                "bitscore": rec["bitscore"],

                "status": status,

            })


            if status in {"HIGH", "CANDIDATE"}:

                accepted.append(
                    (
                        gene,
                        rec["target"],
                        start,
                        end,
                    )
                )


        else:

            key_rows.append({

                "genome": genome,

                "gene": gene,

                "target": "NA",

                "identity": 0,

                "qcov": 0,

                "start": 0,

                "end": 0,

                "strand": "NA",

                "evalue": 1,

                "bitscore": 0,

                "status": "NO_HIT",

            })


    # --------------------------------------------------------
    # Architecture
    # --------------------------------------------------------

    n_detected = len(accepted)


    all_three = (
        n_detected == 3
    )


    same_contig = False

    span = None

    clustered = False


    if all_three:

        contigs = {
            x[1]
            for x in accepted
        }


        same_contig = (
            len(contigs) == 1
        )


        if same_contig:

            starts = [
                x[2]
                for x in accepted
            ]

            ends = [
                x[3]
                for x in accepted
            ]


            span = (
                max(ends)
                - min(starts)
                + 1
            )


            clustered = (
                span <= 50000
            )


    architecture_rows.append({

        "genome": genome,

        "n_key_detected": n_detected,

        "all_EHI": all_three,

        "same_contig": same_contig,

        "span_bp": (
            span
            if span is not None
            else "NA"
        ),

        "clustered_within_50kb": clustered,

    })



# ============================================================
# WRITE KEY-HIT TABLE
# ============================================================

key_file = results / "VfmEHI_all_genomes.tsv"


with open(key_file, "w") as out:

    out.write(
        "genome\tgene\ttarget_contig\t"
        "identity_pct\tquery_coverage_pct\t"
        "start\tend\tstrand\t"
        "evalue\tbitscore\tclassification\n"
    )


    for r in key_rows:

        out.write(
            f"{r['genome']}\t"
            f"{r['gene']}\t"
            f"{r['target']}\t"
            f"{r['identity']:.2f}\t"
            f"{r['qcov']:.1f}\t"
            f"{r['start']}\t"
            f"{r['end']}\t"
            f"{r['strand']}\t"
            f"{r['evalue']:.3g}\t"
            f"{r['bitscore']:.1f}\t"
            f"{r['status']}\n"
        )



# ============================================================
# WRITE ARCHITECTURE TABLE
# ============================================================

architecture_file = (
    results / "VfmEHI_architecture_summary.tsv"
)


with open(architecture_file, "w") as out:

    out.write(
        "genome\t"
        "VfmEHI_detected\t"
        "all_three_EHI\t"
        "same_contig\t"
        "EHI_span_bp\t"
        "clustered_within_50kb\n"
    )


    for r in architecture_rows:

        out.write(
            f"{r['genome']}\t"
            f"{r['n_key_detected']}\t"
            f"{r['all_EHI']}\t"
            f"{r['same_contig']}\t"
            f"{r['span_bp']}\t"
            f"{r['clustered_within_50kb']}\n"
        )



# ============================================================
# BROADER VFM-REGION SUPPORT
# ============================================================

region_rows = []


for genome in genomes:


    if genome == "Dsolani_MK10":

        path = (
            raw
            / "Vfm_region_vs_Dsolani_MK10.tsv"
        )

    else:

        path = (
            raw
            / f"Vfm_region_vs_{genome}.tsv"
        )


    hits = read_best_hits(path)


    high = 0

    candidate = 0

    weak = 0


    accepted_contigs = defaultdict(list)


    for query, rec in hits.items():


        status = classify(
            rec["identity"],
            rec["qcov"],
            rec["evalue"],
        )


        if status == "HIGH":

            high += 1


        elif status == "CANDIDATE":

            candidate += 1


        else:

            weak += 1


        if status in {
            "HIGH",
            "CANDIDATE",
        }:

            start = min(
                rec["sstart"],
                rec["send"],
            )

            end = max(
                rec["sstart"],
                rec["send"],
            )


            accepted_contigs[
                rec["target"]
            ].append(
                (start, end)
            )


    top_contig = "NA"

    top_count = 0

    top_span = "NA"


    if accepted_contigs:

        top_contig, coords = max(
            accepted_contigs.items(),
            key=lambda x: len(x[1]),
        )


        top_count = len(coords)


        top_span = (
            max(end for start, end in coords)
            -
            min(start for start, end in coords)
            + 1
        )


    region_rows.append({

        "genome": genome,

        "high": high,

        "candidate": candidate,

        "weak": weak,

        "total_detected": high + candidate,

        "top_contig": top_contig,

        "top_contig_count": top_count,

        "top_span": top_span,

    })



region_file = (
    results / "Vfm_region_conservation_summary.tsv"
)


with open(region_file, "w") as out:

    out.write(
        "genome\t"
        "high_confidence_hits\t"
        "candidate_hits\t"
        "weak_hits\t"
        "accepted_region_hits\t"
        "top_contig\t"
        "accepted_hits_on_top_contig\t"
        "top_contig_span_bp\n"
    )


    for r in region_rows:

        out.write(
            f"{r['genome']}\t"
            f"{r['high']}\t"
            f"{r['candidate']}\t"
            f"{r['weak']}\t"
            f"{r['total_detected']}\t"
            f"{r['top_contig']}\t"
            f"{r['top_contig_count']}\t"
            f"{r['top_span']}\n"
        )



# ============================================================
# PRIMARY 20GA0316 TABLE
# ============================================================

primary_file = (
    results / "20GA0316_VfmEHI.tsv"
)


with open(primary_file, "w") as out:

    out.write(
        "gene\ttarget_contig\tidentity_pct\t"
        "query_coverage_pct\tstart\tend\tstrand\t"
        "evalue\tbitscore\tclassification\n"
    )


    for r in key_rows:

        if r["genome"] != primary:

            continue


        out.write(
            f"{r['gene']}\t"
            f"{r['target']}\t"
            f"{r['identity']:.2f}\t"
            f"{r['qcov']:.1f}\t"
            f"{r['start']}\t"
            f"{r['end']}\t"
            f"{r['strand']}\t"
            f"{r['evalue']:.3g}\t"
            f"{r['bitscore']:.1f}\t"
            f"{r['status']}\n"
        )



# ============================================================
# D. SOLANI POSITIVE CONTROL TABLE
# ============================================================

mk10_file = (
    results / "Dsolani_MK10_VfmEHI.tsv"
)


with open(mk10_file, "w") as out:

    out.write(
        "gene\ttarget_contig\tidentity_pct\t"
        "query_coverage_pct\tstart\tend\tstrand\t"
        "evalue\tbitscore\tclassification\n"
    )


    for r in key_rows:

        if r["genome"] != "Dsolani_MK10":

            continue


        out.write(
            f"{r['gene']}\t"
            f"{r['target']}\t"
            f"{r['identity']:.2f}\t"
            f"{r['qcov']:.1f}\t"
            f"{r['start']}\t"
            f"{r['end']}\t"
            f"{r['strand']}\t"
            f"{r['evalue']:.3g}\t"
            f"{r['bitscore']:.1f}\t"
            f"{r['status']}\n"
        )


print()
print("Vfm summary complete.")
print()
print(f"Key hits: {key_file}")
print(f"Architecture: {architecture_file}")
print(f"Region conservation: {region_file}")
print()

PY


    mark_done "08_Vfm_summary"

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
echo " VFM ANALYSIS COMPLETE"
echo "============================================================"
echo

echo "Canonical D. dadantii Vfm key proteins:"
echo "  ${RESULTS}/Ddadantii_Vfm_key_reference.tsv"

echo
echo "Canonical Vfm-region proteins:"
echo "  ${RESULTS}/Ddadantii_Vfm_region_reference.tsv"

echo
echo "D. solani MK10 positive control:"
echo "  ${RESULTS}/Dsolani_MK10_VfmEHI.tsv"

echo
echo "R. badensis 20GA0316:"
echo "  ${RESULTS}/20GA0316_VfmEHI.tsv"

echo
echo "All R. badensis genomes:"
echo "  ${RESULTS}/VfmEHI_all_genomes.tsv"

echo
echo "VfmE/H/I genomic architecture:"
echo "  ${RESULTS}/VfmEHI_architecture_summary.tsv"

echo
echo "Broader Vfm-region conservation:"
echo "  ${RESULTS}/Vfm_region_conservation_summary.tsv"

echo
echo "Checkpoints:"
echo "  ${STATE}"

echo
echo "Finished:"
echo "  $(date)"

echo



# ============================================================
# FUTURE STEPS CAN BE ADDED HERE
# ============================================================
#
# STEP 09:
#   inspect any candidate Rouxiella Vfm proteins against
#   conserved-domain annotations
#
# STEP 10:
#   extract candidate genomic neighborhoods
#
# STEP 11:
#   make D. solani vs R. badensis synteny figure
#
# STEP 12:
#   compare Vfm presence with sol cluster presence
#
# Re-running this same script will skip STEP 01-08.
#
# ============================================================

# ============================================================
# STEP 09
# EXTRACT THE ACTUAL CANONICAL VfmA-VfmZ PROTEINS
# ============================================================
#
# Goal:
#
# Extract only the 26 named Vfm proteins:
#
#   VfmA, VfmB, ... VfmZ
#
# from D. dadantii 3937 (NC_014500.1).
#
# This replaces the earlier broad +/-20-kb neighborhood
# analysis with a gene-specific Vfm presence/absence test.
#
# ============================================================

if ! step_done "09_extract_VfmA_Z"; then

    echo
    echo "============================================================"
    echo "STEP 09: Extracting canonical VfmA-VfmZ proteins"
    echo "============================================================"
    echo


    export DDAD_CDS
    export QUERY
    export RESULTS


    python <<'PY'

from pathlib import Path
import os
import re
import string
import sys


fasta = Path(os.environ["DDAD_CDS"])

query_dir = Path(os.environ["QUERY"])

results_dir = Path(os.environ["RESULTS"])


expected = [
    f"vfm{x}"
    for x in string.ascii_uppercase
]


# ============================================================
# FASTA PARSER
# ============================================================

def read_fasta(path):

    records = []

    header = None
    seq = []


    with open(path) as handle:

        for line in handle:

            line = line.rstrip()


            if line.startswith(">"):

                if header is not None:

                    records.append(
                        (header, "".join(seq))
                    )


                header = line[1:]
                seq = []


            else:

                seq.append(line.strip())


        if header is not None:

            records.append(
                (header, "".join(seq))
            )


    return records



records = read_fasta(fasta)



# ============================================================
# EXTRACT GENE / LOCUS TAG / LOCATION FROM NCBI HEADERS
# ============================================================

def field(header, key):

    m = re.search(
        rf"\[{re.escape(key)}=([^\]]+)\]",
        header,
        flags=re.I,
    )

    return m.group(1) if m else None



def coords(header):

    location = field(header, "location")


    if not location:

        return None


    nums = [
        int(x)
        for x in re.findall(r"\d+", location)
    ]


    if len(nums) < 2:

        return None


    return min(nums), max(nums)



found = {}


for header, seq in records:

    gene = field(header, "gene")


    if gene is None:

        continue


    gene = gene.lower()


    if gene in expected:

        found[gene] = {

            "header": header,

            "sequence": seq,

            "locus_tag": (
                field(header, "locus_tag")
                or "NA"
            ),

            "coords": coords(header),

        }



# ============================================================
# REPORT WHAT NCBI ACTUALLY CONTAINS
# ============================================================

debug = results_dir / "Ddadantii_VfmA_Z_extraction_QC.txt"


with open(debug, "w") as out:

    out.write(
        f"Expected genes: {len(expected)}\n"
    )

    out.write(
        f"Found genes: {len(found)}\n\n"
    )


    for gene in expected:

        if gene in found:

            rec = found[gene]

            out.write(
                f"{gene}\t"
                f"{rec['locus_tag']}\t"
                f"{len(rec['sequence'])} aa\n"
            )

        else:

            out.write(
                f"{gene}\tMISSING\n"
            )



missing = [
    gene
    for gene in expected
    if gene not in found
]


if missing:

    print(
        "ERROR: Could not recover all 26 named Vfm proteins.",
        file=sys.stderr,
    )

    print(
        "Missing:",
        ", ".join(missing),
        file=sys.stderr,
    )

    print(
        f"See: {debug}",
        file=sys.stderr,
    )

    sys.exit(1)



# ============================================================
# DETERMINE TRUE REFERENCE ORDER FROM GENOMIC COORDINATES
#
# Do NOT assume alphabetical order = genomic order.
# ============================================================

ordered = sorted(

    found.items(),

    key=lambda x: (
        x[1]["coords"][0]
        if x[1]["coords"] is not None
        else 10**20
    )

)



# ============================================================
# WRITE 26-PROTEIN QUERY FASTA
# ============================================================

out_fasta = (
    query_dir
    / "Ddadantii_VfmA_Z_26proteins.faa"
)


with open(out_fasta, "w") as out:

    for gene, rec in ordered:

        out.write(
            f">{gene}|{rec['locus_tag']}|Ddadantii3937\n"
        )


        seq = rec["sequence"]


        for i in range(0, len(seq), 70):

            out.write(
                seq[i:i+70] + "\n"
            )



# ============================================================
# WRITE REFERENCE TABLE
# ============================================================

reference_table = (
    results_dir
    / "Ddadantii_VfmA_Z_reference.tsv"
)


with open(reference_table, "w") as out:

    out.write(
        "reference_rank\tgene\tlocus_tag\t"
        "start\tend\tprotein_length\n"
    )


    for rank, (gene, rec) in enumerate(
        ordered,
        start=1,
    ):

        start, end = rec["coords"]


        out.write(
            f"{rank}\t"
            f"{gene}\t"
            f"{rec['locus_tag']}\t"
            f"{start}\t"
            f"{end}\t"
            f"{len(rec['sequence'])}\n"
        )



print()
print("Successfully extracted all 26 Vfm proteins.")
print()

print("Reference genomic order:")

for rank, (gene, rec) in enumerate(
    ordered,
    start=1,
):

    print(
        rank,
        gene,
        rec["locus_tag"],
        rec["coords"],
    )

PY


    [[ -s \
"${QUERY}/Ddadantii_VfmA_Z_26proteins.faa" ]] || \
        die "Canonical 26-protein Vfm FASTA was not produced."


    N_VFM=$(

        grep -c '^>' \
"${QUERY}/Ddadantii_VfmA_Z_26proteins.faa"

    )


    if [[ "${N_VFM}" -ne 26 ]]; then

        die "Expected 26 Vfm proteins but obtained ${N_VFM}."

    fi


    mark_done "09_extract_VfmA_Z"

else

    echo
    echo ">>> STEP 09 already completed — skipping"

fi



# ============================================================
# STEP 10
# SEARCH ALL 26 VFM PROTEINS
# ============================================================
#
# Search:
#
#   D. solani MK10 positive control
#   all R. badensis genomes
#
# ============================================================

if ! step_done "10_search_VfmA_Z"; then

    echo
    echo "============================================================"
    echo "STEP 10: Searching VfmA-VfmZ in all genomes"
    echo "============================================================"
    echo


    VFM26_QUERY="${QUERY}/Ddadantii_VfmA_Z_26proteins.faa"


    # ========================================================
    # D. SOLANI POSITIVE CONTROL
    # ========================================================

    tblastn \
        -query "${VFM26_QUERY}" \
        -db "${DB}/Dsolani_MK10/MK10" \
        -evalue 1e-5 \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -num_threads "${THREADS}" \
        -outfmt \
'6 qseqid sseqid pident length qlen qcovhsp sstart send evalue bitscore' \
        > "${RAW}/VfmA_Z_vs_Dsolani_MK10.tsv"


    # ========================================================
    # ROUXIELLA GENOMES
    # ========================================================

    RB_GENOMES=( "${ROOT}/02_genomes/"*.fna )


    for FNA in "${RB_GENOMES[@]}"; do


        ACC="$(basename "${FNA}" .fna)"


        GENOME_DONE="${STATE}/10_VfmAZ_${ACC}.done"


        if [[ -f "${GENOME_DONE}" ]]; then

            echo
            echo "${ACC}: VfmA-Z search already completed — skipping"

            continue

        fi


        echo
        echo "Searching 26 Vfm proteins in ${ACC}"


        tblastn \
            -query "${VFM26_QUERY}" \
            -db "${DB}/${ACC}/${ACC}" \
            -evalue 1e-5 \
            -max_target_seqs 10 \
            -max_hsps 1 \
            -num_threads "${THREADS}" \
            -outfmt \
'6 qseqid sseqid pident length qlen qcovhsp sstart send evalue bitscore' \
            > "${RAW}/VfmA_Z_vs_${ACC}.tsv"


        touch "${GENOME_DONE}"


    done


    mark_done "10_search_VfmA_Z"

else

    echo
    echo ">>> STEP 10 already completed — skipping"

fi



# ============================================================
# STEP 11
# CREATE VfmA-Z PRESENCE / ABSENCE MATRIX
# AND TEST FOR A COHERENT VFM LOCUS
# ============================================================

if ! step_done "11_VfmA_Z_summary"; then

    echo
    echo "============================================================"
    echo "STEP 11: VfmA-Z presence/absence and synteny analysis"
    echo "============================================================"
    echo


    export ROOT
    export RAW
    export RESULTS
    export RB_PRIMARY


    python <<'PY'

from pathlib import Path
from collections import defaultdict
import os


root = Path(os.environ["ROOT"])

raw = Path(os.environ["RAW"])

results = Path(os.environ["RESULTS"])

primary = os.environ["RB_PRIMARY"]



# ============================================================
# REFERENCE GENE ORDER
# ============================================================

reference_file = (
    results
    / "Ddadantii_VfmA_Z_reference.tsv"
)


reference_order = []

reference_rank = {}


with open(reference_file) as handle:

    next(handle)


    for line in handle:

        fields = line.rstrip().split("\t")


        rank = int(fields[0])

        gene = fields[1]


        reference_order.append(gene)

        reference_rank[gene] = rank



# ============================================================
# GENOMES
# ============================================================

genomes = [
    "Dsolani_MK10"
]


genomes += sorted(

    p.stem
    for p in
    (root / "02_genomes").glob("*.fna")

)



# ============================================================
# CLASSIFICATION
#
# HIGH:
#
#   >= 45% amino-acid identity
#   >= 80% query coverage
#   E <= 1e-20
#
# CANDIDATE:
#
#   >= 30% identity
#   >= 60% query coverage
#   E <= 1e-10
#
# WEAK:
#
#   detectable similarity but below candidate threshold
#
# ABSENT:
#
#   no BLAST hit
#
#
# IMPORTANT:
#
# Even HIGH/CANDIDATE does not by itself prove Vfm orthology.
# Co-localization and synteny are assessed separately.
#
# ============================================================

def classify(identity, qcov, evalue):

    if (
        identity >= 45
        and qcov >= 80
        and evalue <= 1e-20
    ):

        return "HIGH"


    if (
        identity >= 30
        and qcov >= 60
        and evalue <= 1e-10
    ):

        return "CANDIDATE"


    return "WEAK"



# ============================================================
# READ BEST HIT FOR EACH QUERY
# ============================================================

def read_best_hits(path):

    best = {}


    if not path.exists():

        return best


    with open(path) as handle:

        for line in handle:


            if not line.strip():

                continue


            p = line.rstrip().split("\t")


            query_full = p[0]


            gene = query_full.split("|")[0]


            rec = {

                "gene": gene,

                "query_full": query_full,

                "contig": p[1],

                "identity": float(p[2]),

                "alignment_length": int(p[3]),

                "qlen": int(p[4]),

                "qcov": float(p[5]),

                "sstart": int(p[6]),

                "send": int(p[7]),

                "evalue": float(p[8]),

                "bitscore": float(p[9]),

            }


            if (
                gene not in best
                or rec["bitscore"]
                > best[gene]["bitscore"]
            ):

                best[gene] = rec


    return best



# ============================================================
# SYNTENY ORDER SCORE
#
# Compare order of accepted genes on the target contig to
# their true D. dadantii reference order.
#
# Because an entire locus can invert, we evaluate both:
#
#   forward orientation
#   reverse orientation
#
# Score:
#
#   1.0 = perfect order
#   0.5 = partially conserved
#   etc.
#
# ============================================================

def order_concordance(genes):

    if len(genes) < 2:

        return None


    ranks = [
        reference_rank[g]
        for g in genes
    ]


    concordant = 0

    discordant = 0


    for i in range(len(ranks)):

        for j in range(i + 1, len(ranks)):


            if ranks[i] < ranks[j]:

                concordant += 1

            elif ranks[i] > ranks[j]:

                discordant += 1


    total = concordant + discordant


    if total == 0:

        return None


    forward = concordant / total

    reverse = discordant / total


    return max(
        forward,
        reverse,
    )



# ============================================================
# PROCESS ALL GENOMES
# ============================================================

long_rows = []

matrix = {}

cluster_rows = []


for genome in genomes:


    if genome == "Dsolani_MK10":

        path = (
            raw
            / "VfmA_Z_vs_Dsolani_MK10.tsv"
        )

    else:

        path = (
            raw
            / f"VfmA_Z_vs_{genome}.tsv"
        )


    best = read_best_hits(path)


    matrix[genome] = {}


    accepted_by_contig = defaultdict(list)


    n_high = 0

    n_candidate = 0

    n_weak = 0

    n_absent = 0


    for gene in reference_order:


        if gene not in best:


            status = "ABSENT"

            n_absent += 1


            matrix[genome][gene] = status


            long_rows.append({

                "genome": genome,

                "gene": gene,

                "contig": "NA",

                "identity": 0,

                "qcov": 0,

                "start": 0,

                "end": 0,

                "strand": "NA",

                "evalue": 1,

                "bitscore": 0,

                "status": status,

            })


            continue



        rec = best[gene]


        status = classify(
            rec["identity"],
            rec["qcov"],
            rec["evalue"],
        )


        matrix[genome][gene] = status


        if status == "HIGH":

            n_high += 1


        elif status == "CANDIDATE":

            n_candidate += 1


        else:

            n_weak += 1


        start = min(
            rec["sstart"],
            rec["send"],
        )

        end = max(
            rec["sstart"],
            rec["send"],
        )


        strand = (
            "+"
            if rec["sstart"] <= rec["send"]
            else "-"
        )


        long_rows.append({

            "genome": genome,

            "gene": gene,

            "contig": rec["contig"],

            "identity": rec["identity"],

            "qcov": rec["qcov"],

            "start": start,

            "end": end,

            "strand": strand,

            "evalue": rec["evalue"],

            "bitscore": rec["bitscore"],

            "status": status,

        })


        if status in {
            "HIGH",
            "CANDIDATE",
        }:

            accepted_by_contig[
                rec["contig"]
            ].append({

                "gene": gene,

                "start": start,

                "end": end,

            })



    # ========================================================
    # IDENTIFY BEST PUTATIVE VFM CONTIG
    # ========================================================

    top_contig = "NA"

    top_hits = []

    top_span = "NA"

    synteny = "NA"


    if accepted_by_contig:


        top_contig, top_hits = max(

            accepted_by_contig.items(),

            key=lambda x: len(x[1]),

        )


        target_sorted = sorted(

            top_hits,

            key=lambda x: x["start"],

        )


        starts = [
            x["start"]
            for x in target_sorted
        ]


        ends = [
            x["end"]
            for x in target_sorted
        ]


        top_span = (
            max(ends)
            - min(starts)
            + 1
        )


        target_gene_order = [
            x["gene"]
            for x in target_sorted
        ]


        score = order_concordance(
            target_gene_order
        )


        if score is not None:

            synteny = round(
                score,
                3,
            )



    accepted_total = (
        n_high
        + n_candidate
    )


    top_count = len(top_hits)



    # ========================================================
    # COHERENT LOCUS DEFINITION
    #
    # Conservative manuscript-level criterion:
    #
    #   >=20/26 accepted Vfm proteins
    #   >=20 on one contig
    #   locus span <=60 kb
    #   gene-order concordance >=0.80
    #
    # ========================================================

    coherent = False


    if (
        accepted_total >= 20
        and top_count >= 20
        and top_span != "NA"
        and top_span <= 60000
        and synteny != "NA"
        and synteny >= 0.80
    ):

        coherent = True



    cluster_rows.append({

        "genome": genome,

        "high": n_high,

        "candidate": n_candidate,

        "weak": n_weak,

        "absent": n_absent,

        "accepted_total": accepted_total,

        "top_contig": top_contig,

        "top_count": top_count,

        "top_span": top_span,

        "synteny": synteny,

        "coherent": coherent,

    })



# ============================================================
# WRITE LONG-FORM BEST-HIT TABLE
# ============================================================

long_file = (
    results
    / "VfmA_Z_best_hits_all_genomes.tsv"
)


with open(long_file, "w") as out:


    out.write(
        "genome\tgene\ttarget_contig\t"
        "identity_pct\tquery_coverage_pct\t"
        "start\tend\tstrand\t"
        "evalue\tbitscore\tclassification\n"
    )


    for r in long_rows:


        out.write(
            f"{r['genome']}\t"
            f"{r['gene']}\t"
            f"{r['contig']}\t"
            f"{r['identity']:.2f}\t"
            f"{r['qcov']:.1f}\t"
            f"{r['start']}\t"
            f"{r['end']}\t"
            f"{r['strand']}\t"
            f"{r['evalue']:.3g}\t"
            f"{r['bitscore']:.1f}\t"
            f"{r['status']}\n"
        )



# ============================================================
# WRITE PRESENCE / ABSENCE MATRIX
# ============================================================

matrix_file = (
    results
    / "VfmA_Z_presence_matrix.tsv"
)


with open(matrix_file, "w") as out:


    out.write(
        "genome\t"
        + "\t".join(reference_order)
        + "\n"
    )


    for genome in genomes:


        values = [
            matrix[genome][gene]
            for gene in reference_order
        ]


        out.write(
            genome
            + "\t"
            + "\t".join(values)
            + "\n"
        )



# ============================================================
# WRITE CLUSTER/SYNTENY SUMMARY
# ============================================================

cluster_file = (
    results
    / "VfmA_Z_cluster_summary.tsv"
)


with open(cluster_file, "w") as out:


    out.write(
        "genome\t"
        "HIGH\t"
        "CANDIDATE\t"
        "WEAK\t"
        "ABSENT\t"
        "accepted_total\t"
        "top_contig\t"
        "accepted_on_top_contig\t"
        "top_contig_span_bp\t"
        "gene_order_concordance\t"
        "coherent_Vfm_locus\n"
    )


    for r in cluster_rows:


        out.write(
            f"{r['genome']}\t"
            f"{r['high']}\t"
            f"{r['candidate']}\t"
            f"{r['weak']}\t"
            f"{r['absent']}\t"
            f"{r['accepted_total']}\t"
            f"{r['top_contig']}\t"
            f"{r['top_count']}\t"
            f"{r['top_span']}\t"
            f"{r['synteny']}\t"
            f"{r['coherent']}\n"
        )



# ============================================================
# WRITE PRIMARY 20GA0316 TABLE
# ============================================================

primary_file = (
    results
    / "20GA0316_VfmA_Z.tsv"
)


with open(primary_file, "w") as out:


    out.write(
        "gene\ttarget_contig\t"
        "identity_pct\tquery_coverage_pct\t"
        "start\tend\tstrand\t"
        "evalue\tbitscore\tclassification\n"
    )


    for r in long_rows:


        if r["genome"] != primary:

            continue


        out.write(
            f"{r['gene']}\t"
            f"{r['contig']}\t"
            f"{r['identity']:.2f}\t"
            f"{r['qcov']:.1f}\t"
            f"{r['start']}\t"
            f"{r['end']}\t"
            f"{r['strand']}\t"
            f"{r['evalue']:.3g}\t"
            f"{r['bitscore']:.1f}\t"
            f"{r['status']}\n"
        )



# ============================================================
# WRITE D. SOLANI POSITIVE CONTROL TABLE
# ============================================================

mk10_file = (
    results
    / "Dsolani_MK10_VfmA_Z.tsv"
)


with open(mk10_file, "w") as out:


    out.write(
        "gene\ttarget_contig\t"
        "identity_pct\tquery_coverage_pct\t"
        "start\tend\tstrand\t"
        "evalue\tbitscore\tclassification\n"
    )


    for r in long_rows:


        if r["genome"] != "Dsolani_MK10":

            continue


        out.write(
            f"{r['gene']}\t"
            f"{r['contig']}\t"
            f"{r['identity']:.2f}\t"
            f"{r['qcov']:.1f}\t"
            f"{r['start']}\t"
            f"{r['end']}\t"
            f"{r['strand']}\t"
            f"{r['evalue']:.3g}\t"
            f"{r['bitscore']:.1f}\t"
            f"{r['status']}\n"
        )



print()
print("VfmA-Z analysis completed.")
print()
print(f"Presence matrix: {matrix_file}")
print(f"Cluster summary: {cluster_file}")
print(f"20GA0316 table: {primary_file}")
print()

PY


    mark_done "11_VfmA_Z_summary"

else

    echo
    echo ">>> STEP 11 already completed — skipping"

fi



# ============================================================
# STEP 12
# PRINT FINAL MANUSCRIPT-LEVEL SUMMARY
# ============================================================

if ! step_done "12_VfmA_Z_report"; then

    echo
    echo "============================================================"
    echo "STEP 12: Final VfmA-Z report"
    echo "============================================================"
    echo


    echo
    echo "===== CANONICAL D. DADANTII VFM LOCUS ====="

    column -t \
"${RESULTS}/Ddadantii_VfmA_Z_reference.tsv" \
        || cat \
"${RESULTS}/Ddadantii_VfmA_Z_reference.tsv"


    echo
    echo "===== D. SOLANI POSITIVE CONTROL ====="

    column -t \
"${RESULTS}/Dsolani_MK10_VfmA_Z.tsv" \
        || true


    echo
    echo "===== R. BADENSIS 20GA0316 ====="

    column -t \
"${RESULTS}/20GA0316_VfmA_Z.tsv" \
        || true


    echo
    echo "===== VFM LOCUS SUMMARY ====="

    column -t \
"${RESULTS}/VfmA_Z_cluster_summary.tsv" \
        || true


    echo
    echo "===== PRESENCE / ABSENCE MATRIX ====="

    column -t \
"${RESULTS}/VfmA_Z_presence_matrix.tsv" \
        || true


    mark_done "12_VfmA_Z_report"

else

    echo
    echo ">>> STEP 12 already completed — skipping"

fi