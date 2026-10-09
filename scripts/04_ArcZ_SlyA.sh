#!/usr/bin/env bash

#SBATCH --job-name=Roux_ArcZSlyA
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/Roux_ArcZSlyA_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/Roux_ArcZSlyA_%j.err

set -euo pipefail
shopt -s nullglob


# ============================================================
# OBJECTIVE 3
#
# CONSERVATION OF THE ArcZ-SlyA REGULATORY AXIS
#
# Biological question:
#
# Has Rouxiella badensis retained the recently identified
# ArcZ -> SlyA -> sol regulatory pathway described in
# Dickeya solani?
#
# We test separately:
#
#   1. Is ArcZ conserved?
#   2. Is SlyA conserved?
#   3. Is the slyA 5' regulatory region conserved?
#   4. Does the predicted ArcZ-binding site in the D. solani
#      slyA 5' region remain conserved in R. badensis?
#   5. Does ArcZ retain the complementary pairing sequence?
#
# IMPORTANT:
#
# The exact experimentally mutated nucleotide coordinates from
# the recent D. solani study are not hard-coded here.
#
# Instead:
#   - D. solani MK10 is used as the biological reference.
#   - IntaRNA identifies the strongest candidate ArcZ-slyA
#     interaction in MK10.
#   - That MK10 site is then mapped through sequence alignments
#     to all R. badensis genomes.
#
# If the exact experimentally validated coordinates become
# available from the paper/supplement, they can later replace
# the computationally inferred MK10 site.
#
# ============================================================


# ============================================================
# PATHS
# ============================================================

ROOT="/scratch/al98750/Roux"

WORK="${ROOT}/08_ArcZ_SlyA"

REFERENCE="${WORK}/01_reference"

RFAM="${WORK}/02_Rfam_ArcZ"

SLYA="${WORK}/03_SlyA"

ALIGN="${WORK}/04_alignments"

INTERACTION="${WORK}/05_interaction"

RESULTS="${WORK}/06_results"

TMP="${WORK}/tmp"

STATE="${WORK}/.state"


# Curated nonredundant Rouxiella set

SELECTED_ACCESSIONS="${ROOT}/01_NCBI/selected_accessions.txt"

RB_GENOME_DIR="${ROOT}/02_genomes"

RB_BAKTA_DIR="${ROOT}/03_annotation/bakta"


# D. solani MK10

MK10_BAKTA="${ROOT}/05_QS_ExpIR/01_Dsolani_MK10/bakta_annotation"

MK10_FNA="${MK10_BAKTA}/Dsolani_MK10.fna"

MK10_TSV="${MK10_BAKTA}/Dsolani_MK10.tsv"

MK10_FAA="${MK10_BAKTA}/Dsolani_MK10.faa"


# Dedicated software environment

ENV="${ROOT}/envs/arcZ_slyA"


THREADS="${SLURM_CPUS_PER_TASK:-4}"


# ============================================================
# ANALYSIS SETTINGS
# ============================================================

# Rfam family for ArcZ

ARCZ_RFAM="RF00081"


# slyA region:
#
# extract 300 nt before the translation start plus
# first 30 nt of the slyA CDS.
#
# All sequences are written in transcript orientation.

SLYA_UPSTREAM=300

SLYA_CDS_CONTEXT=30


# ============================================================
# CREATE DIRECTORIES
# ============================================================

mkdir -p \
    "${WORK}" \
    "${REFERENCE}" \
    "${RFAM}" \
    "${SLYA}" \
    "${ALIGN}" \
    "${INTERACTION}" \
    "${RESULTS}" \
    "${TMP}" \
    "${STATE}" \
    "${ROOT}/logs"


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
# INITIALIZE CONDA
# ============================================================

source "$(conda info --base)/etc/profile.d/conda.sh"


# ============================================================
# CREATE ANALYSIS ENVIRONMENT IF NECESSARY
# ============================================================

echo
echo "============================================================"
echo " SOFTWARE ENVIRONMENT"
echo "============================================================"
echo


NEED_ENV=0


for PROGRAM in \
    cmsearch \
    mafft \
    IntaRNA \
    makeblastdb \
    blastn

do

    if [[ ! -x "${ENV}/bin/${PROGRAM}" ]]; then

        NEED_ENV=1

    fi

done


if [[ "${NEED_ENV}" -eq 1 ]]; then

    echo "Creating/updating ArcZ-SlyA environment:"
    echo "  ${ENV}"
    echo


    if [[ -d "${ENV}" ]]; then

        conda install \
            -p "${ENV}" \
            -c conda-forge \
            -c bioconda \
            python=3.11 \
            blast \
            infernal \
            mafft \
            intarna \
            -y

    else

        conda create \
            -p "${ENV}" \
            -c conda-forge \
            -c bioconda \
            python=3.11 \
            blast \
            infernal \
            mafft \
            intarna \
            -y

    fi

fi


conda activate "${ENV}"


echo "Using environment:"
echo "  ${CONDA_PREFIX}"

echo

for PROGRAM in \
    python \
    cmsearch \
    mafft \
    IntaRNA \
    blastn

do

    command -v "${PROGRAM}" >/dev/null 2>&1 || \
        die "${PROGRAM} is unavailable."

    echo "FOUND: ${PROGRAM}"

done


# ============================================================
# JOB INFORMATION
# ============================================================

echo
echo "============================================================"
echo " ArcZ-SlyA CONSERVATION ANALYSIS"
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
echo "Working directory:"
echo "  ${WORK}"

echo
echo "Started:"
echo "  $(date)"

echo


# ============================================================
# STEP 01
# INPUT VALIDATION
# ============================================================

if ! step_done "01_inputs"; then

    echo
    echo "============================================================"
    echo "STEP 01: Validating inputs"
    echo "============================================================"
    echo


    [[ -s "${SELECTED_ACCESSIONS}" ]] || \
        die "Missing curated accession list."


    [[ -s "${MK10_FNA}" ]] || \
        die "Missing D. solani MK10 genome."


    [[ -s "${MK10_TSV}" ]] || \
        die "Missing D. solani MK10 Bakta TSV."


    [[ -s "${MK10_FAA}" ]] || \
        die "Missing D. solani MK10 protein FASTA."


    mapfile -t RB_ACCESSIONS < <(

        sed 's/\r$//' "${SELECTED_ACCESSIONS}" |
        sed '/^[[:space:]]*$/d'

    )


    [[ "${#RB_ACCESSIONS[@]}" -gt 0 ]] || \
        die "Curated accession list is empty."


    echo "Curated R. badensis genomes:"
    echo "  ${#RB_ACCESSIONS[@]}"
    echo


    for ACC in "${RB_ACCESSIONS[@]}"; do

        [[ -s "${RB_GENOME_DIR}/${ACC}.fna" ]] || \
            die "Genome missing: ${ACC}"

        [[ -s "${RB_BAKTA_DIR}/${ACC}/${ACC}.tsv" ]] || \
            die "Bakta TSV missing: ${ACC}"

        [[ -s "${RB_BAKTA_DIR}/${ACC}/${ACC}.faa" ]] || \
            die "Bakta protein FASTA missing: ${ACC}"

        echo "FOUND: ${ACC}"

    done


    mark_done "01_inputs"

else

    echo
    echo ">>> STEP 01 already completed — skipping"

fi


# ============================================================
# STEP 02
# DOWNLOAD ArcZ RFAM COVARIANCE MODEL
# ============================================================

if ! step_done "02_ArcZ_Rfam_model"; then

    echo
    echo "============================================================"
    echo "STEP 02: Obtaining Rfam ArcZ model RF00081"
    echo "============================================================"
    echo


    ARCZ_CM="${REFERENCE}/${ARCZ_RFAM}.cm"


    wget -q \
        -O "${ARCZ_CM}" \
        "https://rfam.org/family/${ARCZ_RFAM}/cm"


    [[ -s "${ARCZ_CM}" ]] || \
        die "Rfam ArcZ covariance model download failed."


    grep -q 'INFERNAL' "${ARCZ_CM}" || \
        die "Downloaded Rfam file does not appear to be a CM."


    echo "ArcZ covariance model:"
    echo "  ${ARCZ_CM}"


    mark_done "02_ArcZ_Rfam_model"

else

    echo
    echo ">>> STEP 02 already completed — skipping"

fi


ARCZ_CM="${REFERENCE}/${ARCZ_RFAM}.cm"


# ============================================================
# STEP 03
# SEARCH FOR ArcZ USING RFAM / INFERNAL
# ============================================================

if ! step_done "03_ArcZ_cmsearch"; then

    echo
    echo "============================================================"
    echo "STEP 03: Searching for ArcZ"
    echo "============================================================"
    echo


    # --------------------------------------------------------
    # D. solani MK10
    # --------------------------------------------------------

    MK10_DONE="${STATE}/03_ArcZ_Dsolani_MK10.done"


    if [[ ! -f "${MK10_DONE}" ]]; then

        echo "Searching ArcZ in D. solani MK10"


        cmsearch \
            --cpu "${THREADS}" \
            --cut_ga \
            --tblout "${RFAM}/Dsolani_MK10.tbl" \
            "${ARCZ_CM}" \
            "${MK10_FNA}" \
            > "${RFAM}/Dsolani_MK10.cmsearch.txt"


        touch "${MK10_DONE}"

    fi


    # --------------------------------------------------------
    # R. badensis genomes
    # --------------------------------------------------------

    while IFS= read -r ACC || [[ -n "${ACC}" ]]; do

        ACC="${ACC//$'\r'/}"

        [[ -z "${ACC}" ]] && continue


        DONE="${STATE}/03_ArcZ_${ACC}.done"


        if [[ -f "${DONE}" ]]; then

            echo "${ACC}: ArcZ search done — skipping"

            continue

        fi


        echo "Searching ArcZ in ${ACC}"


        cmsearch \
            --cpu "${THREADS}" \
            --cut_ga \
            --tblout "${RFAM}/${ACC}.tbl" \
            "${ARCZ_CM}" \
            "${RB_GENOME_DIR}/${ACC}.fna" \
            > "${RFAM}/${ACC}.cmsearch.txt"


        touch "${DONE}"


    done < "${SELECTED_ACCESSIONS}"


    mark_done "03_ArcZ_cmsearch"

else

    echo
    echo ">>> STEP 03 already completed — skipping"

fi


# ============================================================
# STEP 04
# EXTRACT BEST ArcZ HIT FROM EACH GENOME
# ============================================================

if ! step_done "04_extract_ArcZ"; then

    echo
    echo "============================================================"
    echo "STEP 04: Extracting ArcZ sequences"
    echo "============================================================"
    echo


    export ROOT
    export RFAM
    export RESULTS
    export SELECTED_ACCESSIONS
    export RB_GENOME_DIR
    export MK10_FNA


    python <<'PY'

from pathlib import Path
import csv
import os


rfam = Path(os.environ["RFAM"])
results = Path(os.environ["RESULTS"])

selected = Path(os.environ["SELECTED_ACCESSIONS"])

rb_genome_dir = Path(os.environ["RB_GENOME_DIR"])

mk10_fna = Path(os.environ["MK10_FNA"])


# ------------------------------------------------------------
# FASTA helpers
# ------------------------------------------------------------

def read_fasta(path):

    records = {}

    name = None
    seq = []

    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    records[name] = "".join(seq).upper()

                name = line[1:].split()[0]
                seq = []

            else:

                seq.append(line)

        if name is not None:
            records[name] = "".join(seq).upper()

    return records


def revcomp(seq):

    table = str.maketrans(
        "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
        "TGCAYRKMSWVHDBNtgcayrkmswvhdbn"
    )

    return seq.translate(table)[::-1]


# ------------------------------------------------------------
# Genome list
# ------------------------------------------------------------

with open(selected) as fh:

    rb_genomes = [
        x.strip()
        for x in fh
        if x.strip()
    ]


genomes = ["Dsolani_MK10"] + rb_genomes


# ------------------------------------------------------------
# Parse Infernal tblout
#
# Expected cmsearch tblout fields include:
#
# target name
# accession
# query name
# accession
# mdl
# mdl from
# mdl to
# seq from
# seq to
# strand
# trunc
# pass
# gc
# bias
# score
# E-value
# inc
#
# ------------------------------------------------------------

summary = []

arc_sequences = {}


for genome in genomes:

    tbl = rfam / f"{genome}.tbl"

    if genome == "Dsolani_MK10":
        fna = mk10_fna
    else:
        fna = rb_genome_dir / f"{genome}.fna"


    candidates = []


    if tbl.exists():

        with open(tbl) as fh:

            for line in fh:

                if not line.strip():
                    continue

                if line.startswith("#"):
                    continue


                p = line.split()


                if len(p) < 17:
                    continue


                try:

                    rec = {

                        "contig": p[0],

                        "seq_from": int(p[7]),

                        "seq_to": int(p[8]),

                        "strand": p[9],

                        "score": float(p[14]),

                        "evalue": float(p[15]),

                    }

                except (ValueError, IndexError):

                    continue


                candidates.append(rec)


    if not candidates:

        summary.append({

            "genome": genome,

            "found": False,

            "contig": "NA",

            "start": 0,

            "end": 0,

            "strand": "NA",

            "length": 0,

            "bitscore": 0,

            "evalue": 1,

        })

        continue


    best = max(
        candidates,
        key=lambda x: x["score"]
    )


    fasta = read_fasta(fna)


    if best["contig"] not in fasta:

        raise RuntimeError(
            f"ArcZ target contig {best['contig']} "
            f"not found for {genome}"
        )


    start = min(
        best["seq_from"],
        best["seq_to"]
    )

    end = max(
        best["seq_from"],
        best["seq_to"]
    )


    seq = fasta[best["contig"]][
        start - 1:end
    ]


    if best["strand"] == "-":
        seq = revcomp(seq)


    arc_sequences[genome] = seq


    summary.append({

        "genome": genome,

        "found": True,

        "contig": best["contig"],

        "start": start,

        "end": end,

        "strand": best["strand"],

        "length": len(seq),

        "bitscore": best["score"],

        "evalue": best["evalue"],

    })


# ------------------------------------------------------------
# Write summary
# ------------------------------------------------------------

summary_file = results / "ArcZ_Rfam_summary.tsv"


with open(summary_file, "w", newline="") as out:

    writer = csv.DictWriter(

        out,

        delimiter="\t",

        fieldnames=[

            "genome",

            "found",

            "contig",

            "start",

            "end",

            "strand",

            "length",

            "bitscore",

            "evalue",

        ]

    )

    writer.writeheader()

    writer.writerows(summary)


# ------------------------------------------------------------
# Write ArcZ FASTA
# ------------------------------------------------------------

fasta_file = results / "ArcZ_all_genomes.fna"


with open(fasta_file, "w") as out:

    for genome in genomes:

        if genome not in arc_sequences:
            continue

        out.write(
            f">{genome}\n"
        )

        seq = arc_sequences[genome]

        for i in range(0, len(seq), 70):

            out.write(
                seq[i:i+70] + "\n"
            )


print()
print("ArcZ extraction complete.")
print(summary_file)
print(fasta_file)
print()

PY


    [[ -s "${RESULTS}/ArcZ_Rfam_summary.tsv" ]] || \
        die "ArcZ summary was not generated."


    [[ -s "${RESULTS}/ArcZ_all_genomes.fna" ]] || \
        die "ArcZ FASTA was not generated."


    mark_done "04_extract_ArcZ"

else

    echo
    echo ">>> STEP 04 already completed — skipping"

fi


# ============================================================
# STEP 05
# IDENTIFY SlyA AND EXTRACT TRANSCRIPT-ORIENTED 5' REGION
# ============================================================

if ! step_done "05_extract_SlyA"; then

    echo
    echo "============================================================"
    echo "STEP 05: Extracting SlyA and slyA 5' regions"
    echo "============================================================"
    echo


    export RESULTS
    export SLYA
    export SELECTED_ACCESSIONS
    export RB_GENOME_DIR
    export RB_BAKTA_DIR
    export MK10_FNA
    export MK10_TSV
    export MK10_FAA
    export SLYA_UPSTREAM
    export SLYA_CDS_CONTEXT


    python <<'PY'

from pathlib import Path
import csv
import os


results = Path(os.environ["RESULTS"])
slya_dir = Path(os.environ["SLYA"])

selected = Path(
    os.environ["SELECTED_ACCESSIONS"]
)

rb_genomes = Path(
    os.environ["RB_GENOME_DIR"]
)

rb_bakta = Path(
    os.environ["RB_BAKTA_DIR"]
)

mk10_fna = Path(
    os.environ["MK10_FNA"]
)

mk10_tsv = Path(
    os.environ["MK10_TSV"]
)

mk10_faa = Path(
    os.environ["MK10_FAA"]
)

upstream_requested = int(
    os.environ["SLYA_UPSTREAM"]
)

cds_context = int(
    os.environ["SLYA_CDS_CONTEXT"]
)


# ------------------------------------------------------------
# FASTA helpers
# ------------------------------------------------------------

def read_fasta(path):

    records = {}

    name = None
    header = None
    seq = []

    with open(path) as fh:

        for line in fh:

            line = line.rstrip()

            if not line:
                continue


            if line.startswith(">"):

                if name is not None:

                    records[name] = {
                        "header": header,
                        "seq": "".join(seq).upper()
                    }


                header = line[1:]
                name = header.split()[0]

                seq = []

            else:

                seq.append(
                    line.strip()
                )


        if name is not None:

            records[name] = {
                "header": header,
                "seq": "".join(seq).upper()
            }


    return records


def revcomp(seq):

    return seq.translate(
        str.maketrans(
            "ACGTRYMKSWBDHVNacgtrymkswbdhvn",
            "TGCAYRKMSWVHDBNtgcayrkmswvhdbn"
        )
    )[::-1]


# ------------------------------------------------------------
# Genome definitions
# ------------------------------------------------------------

with open(selected) as fh:

    rb_accessions = [
        line.strip()
        for line in fh
        if line.strip()
    ]


genomes = ["Dsolani_MK10"] + rb_accessions


summary_rows = []

protein_sequences = {}

window_sequences = {}


# ------------------------------------------------------------
# Process each genome
# ------------------------------------------------------------

for genome in genomes:


    if genome == "Dsolani_MK10":

        tsv = mk10_tsv
        faa = mk10_faa
        fna = mk10_fna

    else:

        tsv = (
            rb_bakta
            / genome
            / f"{genome}.tsv"
        )

        faa = (
            rb_bakta
            / genome
            / f"{genome}.faa"
        )

        fna = (
            rb_genomes
            / f"{genome}.fna"
        )


    slya_hits = []


    with open(tsv) as fh:

        reader = csv.reader(
            fh,
            delimiter="\t"
        )


        for row in reader:

            if len(row) < 8:
                continue


            feature_type = row[1].strip()

            gene = row[6].strip()


            if (
                feature_type.lower() == "cds"
                and gene.lower() == "slya"
            ):

                slya_hits.append(row)


    if not slya_hits:

        summary_rows.append({

            "genome": genome,

            "found": False,

            "contig": "NA",

            "start": 0,

            "end": 0,

            "strand": "NA",

            "locus_tag": "NA",

            "protein_length": 0,

            "upstream_nt": 0,

            "window_length": 0,

        })

        continue


    if len(slya_hits) > 1:

        print(
            f"WARNING: {genome} has "
            f"{len(slya_hits)} slyA annotations; "
            f"using the first."
        )


    row = slya_hits[0]


    contig = row[0]

    start = int(row[2])

    end = int(row[3])

    strand = row[4]

    locus = row[5]


    genome_fasta = read_fasta(fna)


    if contig not in genome_fasta:

        raise RuntimeError(
            f"{genome}: contig {contig} "
            f"not found in genome FASTA."
        )


    contig_seq = genome_fasta[
        contig
    ]["seq"]


    # --------------------------------------------------------
    # Extract SlyA protein
    # --------------------------------------------------------

    protein_fasta = read_fasta(faa)


    protein_record = None


    if locus in protein_fasta:

        protein_record = protein_fasta[
            locus
        ]

    else:

        # fallback: search entire header

        for key, rec in protein_fasta.items():

            if locus in rec["header"]:

                protein_record = rec
                break


    if protein_record is None:

        raise RuntimeError(
            f"{genome}: could not extract "
            f"SlyA protein {locus}"
        )


    protein_seq = protein_record[
        "seq"
    ]


    protein_sequences[
        genome
    ] = protein_seq


    # --------------------------------------------------------
    # Extract:
    #
    #   300 nt upstream
    #   +
    #   first 30 nt of slyA CDS
    #
    # Sequence is oriented as the slyA transcript.
    # --------------------------------------------------------

    if strand == "+":

        translation_start0 = (
            start - 1
        )


        slice_start = max(
            0,
            translation_start0
            - upstream_requested
        )


        slice_end = min(
            len(contig_seq),
            translation_start0
            + cds_context
        )


        window = contig_seq[
            slice_start:slice_end
        ]


        actual_upstream = (
            translation_start0
            - slice_start
        )


    elif strand == "-":

        # The translation start is at genomic "end".
        #
        # Extract CDS context + upstream sequence from
        # the genomic plus strand, then reverse-complement.

        slice_start = max(
            0,
            end - cds_context
        )


        slice_end = min(
            len(contig_seq),
            end + upstream_requested
        )


        window = revcomp(
            contig_seq[
                slice_start:slice_end
            ]
        )


        actual_upstream = (
            slice_end - end
        )


    else:

        raise RuntimeError(
            f"{genome}: invalid strand {strand}"
        )


    window_sequences[
        genome
    ] = window


    summary_rows.append({

        "genome": genome,

        "found": True,

        "contig": contig,

        "start": start,

        "end": end,

        "strand": strand,

        "locus_tag": locus,

        "protein_length": len(
            protein_seq
        ),

        "upstream_nt": actual_upstream,

        "window_length": len(
            window
        ),

    })


# ------------------------------------------------------------
# Write annotation summary
# ------------------------------------------------------------

summary_file = (
    results
    / "SlyA_annotation_summary.tsv"
)


with open(
    summary_file,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(

        out,

        delimiter="\t",

        fieldnames=[

            "genome",

            "found",

            "contig",

            "start",

            "end",

            "strand",

            "locus_tag",

            "protein_length",

            "upstream_nt",

            "window_length",

        ]

    )

    writer.writeheader()

    writer.writerows(
        summary_rows
    )


# ------------------------------------------------------------
# Write SlyA protein FASTA
# ------------------------------------------------------------

protein_file = (
    slya_dir
    / "SlyA_all_genomes.faa"
)


with open(protein_file, "w") as out:

    for genome in genomes:

        if genome not in protein_sequences:
            continue


        out.write(
            f">{genome}\n"
        )


        seq = protein_sequences[
            genome
        ]


        for i in range(
            0,
            len(seq),
            70
        ):

            out.write(
                seq[i:i+70]
                + "\n"
            )


# ------------------------------------------------------------
# Write transcript-oriented slyA 5' windows
# ------------------------------------------------------------

window_file = (
    slya_dir
    / "SlyA_5prime_windows.fna"
)


with open(window_file, "w") as out:

    for genome in genomes:

        if genome not in window_sequences:
            continue


        out.write(
            f">{genome}\n"
        )


        seq = window_sequences[
            genome
        ]


        for i in range(
            0,
            len(seq),
            70
        ):

            out.write(
                seq[i:i+70]
                + "\n"
            )


print()
print("SlyA extraction complete.")
print(summary_file)
print(protein_file)
print(window_file)
print()

PY


    [[ -s "${RESULTS}/SlyA_annotation_summary.tsv" ]] || \
        die "SlyA annotation summary missing."


    [[ -s "${SLYA}/SlyA_all_genomes.faa" ]] || \
        die "SlyA proteins were not extracted."


    [[ -s "${SLYA}/SlyA_5prime_windows.fna" ]] || \
        die "slyA 5' windows were not extracted."


    mark_done "05_extract_SlyA"

else

    echo
    echo ">>> STEP 05 already completed — skipping"

fi


# ============================================================
# STEP 06
# ALIGN ArcZ, SlyA PROTEIN, AND slyA 5' REGION
# ============================================================

if ! step_done "06_alignments"; then

    echo
    echo "============================================================"
    echo "STEP 06: Multiple sequence alignments"
    echo "============================================================"
    echo


    mafft \
        --auto \
        --thread "${THREADS}" \
        "${RESULTS}/ArcZ_all_genomes.fna" \
        > "${ALIGN}/ArcZ.mafft.fna"


    mafft \
        --auto \
        --thread "${THREADS}" \
        "${SLYA}/SlyA_all_genomes.faa" \
        > "${ALIGN}/SlyA_protein.mafft.faa"


    mafft \
        --auto \
        --thread "${THREADS}" \
        "${SLYA}/SlyA_5prime_windows.fna" \
        > "${ALIGN}/SlyA_5prime.mafft.fna"


    mark_done "06_alignments"

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi


# ============================================================
# STEP 07
# CALCULATE FULL-SEQUENCE CONSERVATION AGAINST D. SOLANI
# ============================================================

if ! step_done "07_sequence_conservation"; then

    echo
    echo "============================================================"
    echo "STEP 07: Sequence conservation relative to D. solani"
    echo "============================================================"
    echo


    export ALIGN
    export RESULTS
    export SELECTED_ACCESSIONS


    python <<'PY'

from pathlib import Path
import csv
import os


align_dir = Path(
    os.environ["ALIGN"]
)

results = Path(
    os.environ["RESULTS"]
)

selected = Path(
    os.environ["SELECTED_ACCESSIONS"]
)


REFERENCE = "Dsolani_MK10"


def read_fasta(path):

    records = {}

    name = None
    seq = []


    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if name is not None:

                    records[name] = (
                        "".join(seq)
                        .upper()
                    )


                name = (
                    line[1:]
                    .split()[0]
                )

                seq = []

            else:

                seq.append(line)


        if name is not None:

            records[name] = (
                "".join(seq)
                .upper()
            )


    return records


def identity_to_reference(
    ref,
    query
):

    if len(ref) != len(query):

        raise ValueError(
            "Alignment lengths differ."
        )


    ref_positions = 0

    query_present = 0

    matches = 0


    for r, q in zip(
        ref,
        query
    ):

        if r == "-":
            continue


        ref_positions += 1


        if q != "-":

            query_present += 1


            if r == q:

                matches += 1


    if ref_positions == 0:

        return None, None


    identity = (
        100.0
        * matches
        / ref_positions
    )


    coverage = (
        100.0
        * query_present
        / ref_positions
    )


    return identity, coverage


with open(selected) as fh:

    genomes = [
        "Dsolani_MK10"
    ] + [
        x.strip()
        for x in fh
        if x.strip()
    ]


datasets = {

    "ArcZ":
        align_dir
        / "ArcZ.mafft.fna",

    "SlyA_protein":
        align_dir
        / "SlyA_protein.mafft.faa",

    "slyA_5prime":
        align_dir
        / "SlyA_5prime.mafft.fna",

}


aligned = {

    key: read_fasta(path)

    for key, path in datasets.items()

}


rows = []


for genome in genomes:

    row = {
        "genome": genome
    }


    for key in (
        "ArcZ",
        "SlyA_protein",
        "slyA_5prime",
    ):

        records = aligned[key]


        if (
            REFERENCE not in records
            or genome not in records
        ):

            row[
                f"{key}_identity_pct"
            ] = "NA"

            row[
                f"{key}_coverage_pct"
            ] = "NA"

            continue


        ident, cov = identity_to_reference(

            records[REFERENCE],

            records[genome]

        )


        row[
            f"{key}_identity_pct"
        ] = round(
            ident,
            2
        )


        row[
            f"{key}_coverage_pct"
        ] = round(
            cov,
            2
        )


    rows.append(row)


outfile = (
    results
    / "sequence_conservation_vs_Dsolani.tsv"
)


fields = [

    "genome",

    "ArcZ_identity_pct",

    "ArcZ_coverage_pct",

    "SlyA_protein_identity_pct",

    "SlyA_protein_coverage_pct",

    "slyA_5prime_identity_pct",

    "slyA_5prime_coverage_pct",

]


with open(
    outfile,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(

        out,

        delimiter="\t",

        fieldnames=fields

    )

    writer.writeheader()

    writer.writerows(rows)


print()
print(outfile)
print()

PY


    mark_done "07_sequence_conservation"

else

    echo
    echo ">>> STEP 07 already completed — skipping"

fi


# ============================================================
# STEP 08
# PREPARE SINGLE-GENOME FASTA FILES FOR IntaRNA
# ============================================================

if ! step_done "08_prepare_IntaRNA"; then

    echo
    echo "============================================================"
    echo "STEP 08: Preparing ArcZ-slyA interaction inputs"
    echo "============================================================"
    echo


    mkdir -p \
        "${TMP}/ArcZ_single" \
        "${TMP}/SlyA_single"


    export RESULTS
    export SLYA
    export TMP


    python <<'PY'

from pathlib import Path
import os


results = Path(
    os.environ["RESULTS"]
)

slya = Path(
    os.environ["SLYA"]
)

tmp = Path(
    os.environ["TMP"]
)


def read_fasta(path):

    records = {}

    name = None
    seq = []


    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if name is not None:
                    records[name] = "".join(seq)


                name = (
                    line[1:]
                    .split()[0]
                )

                seq = []

            else:

                seq.append(line)


        if name is not None:
            records[name] = "".join(seq)


    return records


arc = read_fasta(
    results
    / "ArcZ_all_genomes.fna"
)

sly = read_fasta(
    slya
    / "SlyA_5prime_windows.fna"
)


for genome, seq in arc.items():

    path = (
        tmp
        / "ArcZ_single"
        / f"{genome}.fna"
    )


    with open(path, "w") as out:

        out.write(
            f">{genome}_ArcZ\n"
            f"{seq}\n"
        )


for genome, seq in sly.items():

    path = (
        tmp
        / "SlyA_single"
        / f"{genome}.fna"
    )


    with open(path, "w") as out:

        out.write(
            f">{genome}_slyA5prime\n"
            f"{seq}\n"
        )


print("IntaRNA input files prepared.")

PY


    mark_done "08_prepare_IntaRNA"

else

    echo
    echo ">>> STEP 08 already completed — skipping"

fi


# ============================================================
# STEP 09
# PREDICT ArcZ-slyA RNA-RNA INTERACTIONS
# ============================================================

if ! step_done "09_IntaRNA"; then

    echo
    echo "============================================================"
    echo "STEP 09: ArcZ-slyA interaction prediction"
    echo "============================================================"
    echo


    mapfile -t ALL_GENOMES < <(

        {
            echo "Dsolani_MK10"
            sed 's/\r$//' "${SELECTED_ACCESSIONS}" |
                sed '/^[[:space:]]*$/d'
        }

    )


    for GENOME in "${ALL_GENOMES[@]}"; do


        ARC_FILE="${TMP}/ArcZ_single/${GENOME}.fna"

        SLY_FILE="${TMP}/SlyA_single/${GENOME}.fna"


        if [[ ! -s "${ARC_FILE}" ]]; then

            echo "${GENOME}: no ArcZ sequence — skipping IntaRNA"

            continue

        fi


        if [[ ! -s "${SLY_FILE}" ]]; then

            echo "${GENOME}: no slyA 5' sequence — skipping IntaRNA"

            continue

        fi


        DONE="${STATE}/09_IntaRNA_${GENOME}.done"


        if [[ -f "${DONE}" ]]; then

            echo "${GENOME}: IntaRNA already done — skipping"

            continue

        fi


        echo
        echo "Predicting ArcZ-slyA pairing:"
        echo "  ${GENOME}"


        IntaRNA \
            --personality=IntaRNAsTar \
            -t "${SLY_FILE}" \
            -q "${ARC_FILE}" \
            -n 1 \
            --outMode=C \
            --outCsvCols=id1,id2,start1,end1,start2,end2,subseq1,subseq2,hybridDP,E \
            --out="${INTERACTION}/${GENOME}.csv" \
            --default-log-file="${INTERACTION}/${GENOME}.log"


        touch "${DONE}"


    done


    mark_done "09_IntaRNA"

else

    echo
    echo ">>> STEP 09 already completed — skipping"

fi


# ============================================================
# STEP 10
# MAP THE D. SOLANI INTERACTION SITE INTO ROUXIELLA ALIGNMENTS
# ============================================================

if ! step_done "10_interaction_site_conservation"; then

    echo
    echo "============================================================"
    echo "STEP 10: Comparing ArcZ-slyA pairing sites"
    echo "============================================================"
    echo


    export ALIGN
    export INTERACTION
    export RESULTS
    export SELECTED_ACCESSIONS


    python <<'PY'

from pathlib import Path
import csv
import os


align_dir = Path(
    os.environ["ALIGN"]
)

interaction_dir = Path(
    os.environ["INTERACTION"]
)

results = Path(
    os.environ["RESULTS"]
)

selected = Path(
    os.environ["SELECTED_ACCESSIONS"]
)


REFERENCE = "Dsolani_MK10"


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def read_fasta(path):

    records = {}

    name = None
    seq = []


    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue


            if line.startswith(">"):

                if name is not None:
                    records[name] = "".join(seq).upper()


                name = line[1:].split()[0]
                seq = []

            else:

                seq.append(line)


        if name is not None:
            records[name] = "".join(seq).upper()


    return records


def read_intarna(path):

    if not path.exists():
        return None


    with open(path) as fh:

        reader = csv.DictReader(
            fh,
            delimiter=";"
        )


        for row in reader:

            return row


    return None


def alignment_columns_for_ungapped_range(
    aligned_reference,
    start,
    end
):

    position = 0

    columns = []


    for col, nt in enumerate(
        aligned_reference
    ):

        if nt != "-":

            position += 1


        if (
            nt != "-"
            and start <= position <= end
        ):

            columns.append(col)


    return columns


def site_sequence(
    aligned_seq,
    columns
):

    return "".join(
        aligned_seq[i]
        for i in columns
    )


def site_identity(
    ref_seq,
    query_seq
):

    if len(ref_seq) == 0:
        return None


    matches = sum(

        1

        for r, q in zip(
            ref_seq,
            query_seq
        )

        if r == q

    )


    return (
        100.0
        * matches
        / len(ref_seq)
    )


# ------------------------------------------------------------
# Genome list
# ------------------------------------------------------------

with open(selected) as fh:

    genomes = [
        REFERENCE
    ] + [
        x.strip()
        for x in fh
        if x.strip()
    ]


# ------------------------------------------------------------
# Read alignments
# ------------------------------------------------------------

arc_aln = read_fasta(

    align_dir
    / "ArcZ.mafft.fna"

)

sly_aln = read_fasta(

    align_dir
    / "SlyA_5prime.mafft.fna"

)


if REFERENCE not in arc_aln:

    raise RuntimeError(
        "D. solani ArcZ missing from alignment."
    )


if REFERENCE not in sly_aln:

    raise RuntimeError(
        "D. solani slyA 5' region missing from alignment."
    )


# ------------------------------------------------------------
# Read D. solani predicted interaction
# ------------------------------------------------------------

ref_prediction = read_intarna(

    interaction_dir
    / f"{REFERENCE}.csv"

)


if ref_prediction is None:

    raise RuntimeError(
        "No D. solani MK10 IntaRNA interaction was produced."
    )


target_start = int(
    ref_prediction["start1"]
)

target_end = int(
    ref_prediction["end1"]
)

arc_start = int(
    ref_prediction["start2"]
)

arc_end = int(
    ref_prediction["end2"]
)


# ------------------------------------------------------------
# Map D. solani interaction positions to alignment columns
# ------------------------------------------------------------

target_columns = (
    alignment_columns_for_ungapped_range(

        sly_aln[REFERENCE],

        target_start,

        target_end

    )
)


arc_columns = (
    alignment_columns_for_ungapped_range(

        arc_aln[REFERENCE],

        arc_start,

        arc_end

    )
)


ref_target_site = site_sequence(

    sly_aln[REFERENCE],

    target_columns

)


ref_arc_site = site_sequence(

    arc_aln[REFERENCE],

    arc_columns

)


# ------------------------------------------------------------
# Write reference interaction details
# ------------------------------------------------------------

reference_out = (

    results
    / "Dsolani_MK10_predicted_ArcZ_slyA_interaction.tsv"

)


with open(
    reference_out,
    "w"
) as out:

    out.write(
        "parameter\tvalue\n"
    )

    out.write(
        f"slyA_window_start\t{target_start}\n"
    )

    out.write(
        f"slyA_window_end\t{target_end}\n"
    )

    out.write(
        f"ArcZ_start\t{arc_start}\n"
    )

    out.write(
        f"ArcZ_end\t{arc_end}\n"
    )

    out.write(
        f"interaction_energy_kcal_mol\t"
        f"{ref_prediction['E']}\n"
    )

    out.write(
        f"slyA_interaction_sequence\t"
        f"{ref_prediction['subseq1']}\n"
    )

    out.write(
        f"ArcZ_interaction_sequence\t"
        f"{ref_prediction['subseq2']}\n"
    )

    out.write(
        f"hybrid_structure\t"
        f"{ref_prediction['hybridDP']}\n"
    )


# ------------------------------------------------------------
# Compare homologous site across genomes
# ------------------------------------------------------------

rows = []


for genome in genomes:

    prediction = read_intarna(

        interaction_dir
        / f"{genome}.csv"

    )


    row = {

        "genome": genome,

        "MK10_slyA_site_identity_pct": "NA",

        "MK10_ArcZ_site_identity_pct": "NA",

        "MK10_slyA_site_aligned_sequence": "NA",

        "MK10_ArcZ_site_aligned_sequence": "NA",

        "IntaRNA_target_start": "NA",

        "IntaRNA_target_end": "NA",

        "IntaRNA_ArcZ_start": "NA",

        "IntaRNA_ArcZ_end": "NA",

        "IntaRNA_energy_kcal_mol": "NA",

    }


    if genome in sly_aln:

        seq = site_sequence(

            sly_aln[genome],

            target_columns

        )


        row[
            "MK10_slyA_site_aligned_sequence"
        ] = seq


        row[
            "MK10_slyA_site_identity_pct"
        ] = round(

            site_identity(
                ref_target_site,
                seq
            ),

            2

        )


    if genome in arc_aln:

        seq = site_sequence(

            arc_aln[genome],

            arc_columns

        )


        row[
            "MK10_ArcZ_site_aligned_sequence"
        ] = seq


        row[
            "MK10_ArcZ_site_identity_pct"
        ] = round(

            site_identity(
                ref_arc_site,
                seq
            ),

            2

        )


    if prediction is not None:

        row[
            "IntaRNA_target_start"
        ] = prediction["start1"]

        row[
            "IntaRNA_target_end"
        ] = prediction["end1"]

        row[
            "IntaRNA_ArcZ_start"
        ] = prediction["start2"]

        row[
            "IntaRNA_ArcZ_end"
        ] = prediction["end2"]

        row[
            "IntaRNA_energy_kcal_mol"
        ] = prediction["E"]


    rows.append(row)


outfile = (

    results
    / "ArcZ_slyA_interaction_site_conservation.tsv"

)


fields = [

    "genome",

    "MK10_slyA_site_identity_pct",

    "MK10_ArcZ_site_identity_pct",

    "MK10_slyA_site_aligned_sequence",

    "MK10_ArcZ_site_aligned_sequence",

    "IntaRNA_target_start",

    "IntaRNA_target_end",

    "IntaRNA_ArcZ_start",

    "IntaRNA_ArcZ_end",

    "IntaRNA_energy_kcal_mol",

]


with open(
    outfile,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(

        out,

        delimiter="\t",

        fieldnames=fields

    )

    writer.writeheader()

    writer.writerows(rows)


print()
print("Reference predicted interaction:")
print(reference_out)

print()
print("Interaction-site conservation:")
print(outfile)
print()

PY


    mark_done "10_interaction_site_conservation"

else

    echo
    echo ">>> STEP 10 already completed — skipping"

fi


# ============================================================
# STEP 11
# CREATE MASTER CONSERVATION TABLE
# ============================================================

if ! step_done "11_master_summary"; then

    echo
    echo "============================================================"
    echo "STEP 11: Creating final ArcZ-SlyA summary"
    echo "============================================================"
    echo


    export RESULTS
    export SELECTED_ACCESSIONS


    python <<'PY'

from pathlib import Path
import csv
import os


results = Path(
    os.environ["RESULTS"]
)

selected = Path(
    os.environ["SELECTED_ACCESSIONS"]
)


def read_table(
    path,
    key="genome"
):

    data = {}


    with open(path) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )


        for row in reader:

            data[
                row[key]
            ] = row


    return data


arc = read_table(

    results
    / "ArcZ_Rfam_summary.tsv"

)

slya = read_table(

    results
    / "SlyA_annotation_summary.tsv"

)

seq = read_table(

    results
    / "sequence_conservation_vs_Dsolani.tsv"

)

site = read_table(

    results
    / "ArcZ_slyA_interaction_site_conservation.tsv"

)


with open(selected) as fh:

    genomes = [
        "Dsolani_MK10"
    ] + [
        x.strip()
        for x in fh
        if x.strip()
    ]


rows = []


for genome in genomes:

    a = arc.get(
        genome,
        {}
    )

    s = slya.get(
        genome,
        {}
    )

    q = seq.get(
        genome,
        {}
    )

    i = site.get(
        genome,
        {}
    )


    rows.append({

        "genome":
            genome,

        "ArcZ_detected":
            a.get(
                "found",
                "NA"
            ),

        "ArcZ_Rfam_bitscore":
            a.get(
                "bitscore",
                "NA"
            ),

        "ArcZ_Rfam_evalue":
            a.get(
                "evalue",
                "NA"
            ),

        "ArcZ_identity_to_MK10_pct":
            q.get(
                "ArcZ_identity_pct",
                "NA"
            ),

        "SlyA_detected":
            s.get(
                "found",
                "NA"
            ),

        "SlyA_protein_identity_to_MK10_pct":
            q.get(
                "SlyA_protein_identity_pct",
                "NA"
            ),

        "slyA_5prime_identity_to_MK10_pct":
            q.get(
                "slyA_5prime_identity_pct",
                "NA"
            ),

        "MK10_ArcZ_pairing_site_identity_pct":
            i.get(
                "MK10_ArcZ_site_identity_pct",
                "NA"
            ),

        "MK10_slyA_pairing_site_identity_pct":
            i.get(
                "MK10_slyA_site_identity_pct",
                "NA"
            ),

        "IntaRNA_energy_kcal_mol":
            i.get(
                "IntaRNA_energy_kcal_mol",
                "NA"
            ),

        "IntaRNA_target_start":
            i.get(
                "IntaRNA_target_start",
                "NA"
            ),

        "IntaRNA_target_end":
            i.get(
                "IntaRNA_target_end",
                "NA"
            ),

        "IntaRNA_ArcZ_start":
            i.get(
                "IntaRNA_ArcZ_start",
                "NA"
            ),

        "IntaRNA_ArcZ_end":
            i.get(
                "IntaRNA_ArcZ_end",
                "NA"
            ),

    })


outfile = (

    results
    / "ArcZ_SlyA_conservation_summary.tsv"

)


fields = list(
    rows[0].keys()
)


with open(
    outfile,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(

        out,

        delimiter="\t",

        fieldnames=fields

    )

    writer.writeheader()

    writer.writerows(rows)


print()
print(outfile)
print()

PY


    mark_done "11_master_summary"

else

    echo
    echo ">>> STEP 11 already completed — skipping"

fi


# ============================================================
# FINAL REPORT
# ============================================================

echo
echo
echo "============================================================"
echo " ArcZ-SlyA ANALYSIS COMPLETE"
echo "============================================================"
echo

echo "ArcZ Rfam results:"
echo "  ${RESULTS}/ArcZ_Rfam_summary.tsv"

echo
echo "SlyA annotations:"
echo "  ${RESULTS}/SlyA_annotation_summary.tsv"

echo
echo "Full-sequence conservation:"
echo "  ${RESULTS}/sequence_conservation_vs_Dsolani.tsv"

echo
echo "D. solani predicted ArcZ-slyA site:"
echo "  ${RESULTS}/Dsolani_MK10_predicted_ArcZ_slyA_interaction.tsv"

echo
echo "Interaction-site conservation:"
echo "  ${RESULTS}/ArcZ_slyA_interaction_site_conservation.tsv"

echo
echo "MASTER SUMMARY:"
echo "  ${RESULTS}/ArcZ_SlyA_conservation_summary.tsv"

echo
echo "Alignments:"
echo "  ${ALIGN}/ArcZ.mafft.fna"
echo "  ${ALIGN}/SlyA_protein.mafft.faa"
echo "  ${ALIGN}/SlyA_5prime.mafft.fna"

echo
echo "Finished:"
echo "  $(date)"

echo