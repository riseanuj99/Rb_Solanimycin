#!/bin/bash
#SBATCH --job-name=SlyA_sol
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/SlyA_sol_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/SlyA_sol_%j.err

set -euo pipefail


# ============================================================
# OBJECTIVE:
#
# Test conservation of the experimentally demonstrated
# SlyA -> Psol regulatory connection from Dickeya solani
# D s0432-1 in Rouxiella badensis.
#
# Experimental reference:
#   Dickeya solani D s0432-1
#   chromosome: CP017453.1
#
# Psol EMSA primers from Table S5:
#
#   L2256-FAM
#   GAATTGTTTGGCGCGGATATG
#
#   L2257
#   CACAATTTACGATATAACCAGAATGAAAAATC
#
# IMPORTANT:
# Every step has its own checkpoint under .state/.
# Completed steps WILL NOT rerun.
#
# To rerun one step intentionally:
#   rm .state/XX_name.done
#
# ============================================================


# ============================================================
# PATHS
# ============================================================

ROOT="/scratch/al98750/Roux/09_SlyA_sol_exact"

REF="$ROOT/01_reference"
PRIMER="$ROOT/02_primer_mapping"
PROBE="$ROOT/03_Psol_EMSA_probe"
RB="$ROOT/04_Rouxiella_comparison"
RESULTS="$ROOT/06_results"
TMP="$ROOT/tmp"
STATE="$ROOT/.state"

LOGS="/scratch/al98750/Roux/logs"

ENV="/scratch/al98750/Roux/envs/arcZ_slyA"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

ACCESSION="CP017453.1"

GENOME="$REF/Dsolani_Ds0432_1_CP017453.1.fna"


mkdir -p \
    "$REF" \
    "$PRIMER" \
    "$PROBE" \
    "$RB" \
    "$RESULTS" \
    "$TMP" \
    "$STATE" \
    "$LOGS"


# ============================================================
# ENVIRONMENT
# ============================================================

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV"


echo
echo "============================================================"
echo " SlyA -> sol regulatory conservation analysis"
echo "============================================================"
echo
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Node: $(hostname)"
echo "Threads: $THREADS"
echo "Root: $ROOT"
echo "Started: $(date)"
echo


for X in python makeblastdb blastn
do
    command -v "$X" >/dev/null 2>&1 || {
        echo "ERROR: required program not found: $X"
        exit 1
    }

    echo "FOUND: $X"
done


# ============================================================
# STEP 01
# GET D. SOLANI D s0432-1 REFERENCE GENOME
# ============================================================

if [[ ! -f "$STATE/01_reference_genome.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 01: Obtaining D. solani D s0432-1 genome"
    echo "============================================================"
    echo


    if [[ -s "$GENOME" ]]
    then

        echo "Reference genome already exists:"
        echo "$GENOME"

    else

        echo "Downloading $ACCESSION from NCBI..."

        URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=fasta&retmode=text"


        if command -v curl >/dev/null 2>&1
        then

            curl \
                -L \
                --fail \
                --retry 3 \
                "$URL" \
                -o "$GENOME"

        elif command -v wget >/dev/null 2>&1
        then

            wget \
                -O "$GENOME" \
                "$URL"

        else

            echo "ERROR:"
            echo "Neither curl nor wget is available."
            exit 1

        fi

    fi


    # --------------------------------------------------------
    # Basic QC
    # --------------------------------------------------------

    if [[ ! -s "$GENOME" ]]
    then
        echo "ERROR: genome download failed."
        exit 1
    fi


    if ! grep -q '^>' "$GENOME"
    then
        echo "ERROR: downloaded file is not FASTA."
        exit 1
    fi


    GENOME_LENGTH=$(
        grep -v '^>' "$GENOME" \
        | tr -d '\n\r ' \
        | wc -c
    )


    echo
    echo "Genome length:"
    echo "$GENOME_LENGTH bp"


    if [[ "$GENOME_LENGTH" -lt 4000000 ]]
    then

        echo "ERROR:"
        echo "Genome is unexpectedly short."
        echo "Expected approximately 4.9 Mb."

        exit 1
    fi


    echo
    echo "Genome FASTA:"
    echo "$GENOME"


    date > "$STATE/01_reference_genome.done"

else

    echo
    echo ">>> STEP 01 already completed — skipping"

fi


# ============================================================
# STEP 02
# CREATE EXACT EMSA PRIMER FILE
# ============================================================

if [[ ! -f "$STATE/02_EMSA_primers.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 02: Preparing experimental Psol EMSA primers"
    echo "============================================================"
    echo


    cat > "$PRIMER/Psol_EMSA_primers.fna" <<'EOF'
>L2256_FAM
GAATTGTTTGGCGCGGATATG
>L2257
CACAATTTACGATATAACCAGAATGAAAAATC
EOF


    echo "Primer sequences:"
    cat "$PRIMER/Psol_EMSA_primers.fna"


    date > "$STATE/02_EMSA_primers.done"

else

    echo
    echo ">>> STEP 02 already completed — skipping"

fi


# ============================================================
# STEP 03
# MAP EMSA PRIMERS TO D. SOLANI D s0432-1
# ============================================================

if [[ ! -f "$STATE/03_map_EMSA_primers.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 03: Mapping Psol EMSA primers"
    echo "============================================================"
    echo


    DB="$TMP/Ds0432_db"

    makeblastdb \
        -in "$GENOME" \
        -dbtype nucl \
        -out "$DB" \
        > "$TMP/makeblastdb.log"


    blastn \
        -task blastn-short \
        -query "$PRIMER/Psol_EMSA_primers.fna" \
        -db "$DB" \
        -strand both \
        -evalue 1000 \
        -word_size 7 \
        -dust no \
        -num_threads "$THREADS" \
        -outfmt \
'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qseq sseq' \
        > "$PRIMER/Psol_EMSA_primer_hits.tsv"


    echo
    echo "Primer BLAST hits:"
    echo

    column -t -s $'\t' \
        "$PRIMER/Psol_EMSA_primer_hits.tsv"


    if [[ ! -s "$PRIMER/Psol_EMSA_primer_hits.tsv" ]]
    then

        echo "ERROR:"
        echo "No primer matches found."

        exit 1
    fi


    date > "$STATE/03_map_EMSA_primers.done"

else

    echo
    echo ">>> STEP 03 already completed — skipping"

fi


# ============================================================
# STEP 04
# IDENTIFY THE PCR PRIMER PAIR AND EXTRACT EXACT EMSA PROBE
# ============================================================

if [[ ! -f "$STATE/04_extract_EMSA_probe.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 04: Reconstructing exact experimental Psol EMSA probe"
    echo "============================================================"
    echo


    python - \
        "$GENOME" \
        "$PRIMER/Psol_EMSA_primer_hits.tsv" \
        "$PROBE/Psol_EMSA_probe.fna" \
        "$PROBE/Psol_EMSA_probe_coordinates.tsv" <<'PY'

import sys
from pathlib import Path


genome_file = Path(sys.argv[1])
hits_file = Path(sys.argv[2])
fasta_out = Path(sys.argv[3])
coord_out = Path(sys.argv[4])


# ============================================================
# FUNCTIONS
# ============================================================

def read_fasta(path):

    records = {}

    name = None
    seq = []

    for line in path.read_text().splitlines():

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


def rc(seq):

    return seq.translate(
        str.maketrans(
            "ACGTN",
            "TGCAN"
        )
    )[::-1]


# ============================================================
# READ GENOME
# ============================================================

genome = read_fasta(genome_file)


# ============================================================
# READ BLAST HITS
# ============================================================

hits = []


with open(hits_file) as fh:

    for line in fh:

        p = line.rstrip().split("\t")

        if len(p) < 15:
            continue

        hits.append({
            "primer": p[0],
            "contig": p[1],
            "pident": float(p[2]),
            "length": int(p[3]),
            "mismatch": int(p[4]),
            "qstart": int(p[6]),
            "qend": int(p[7]),
            "sstart": int(p[8]),
            "send": int(p[9]),
            "strand": p[10],
            "bitscore": float(p[12])
        })


p1 = [
    h for h in hits
    if h["primer"] == "L2256_FAM"
]

p2 = [
    h for h in hits
    if h["primer"] == "L2257"
]


if not p1 or not p2:

    raise RuntimeError(
        "Could not identify both EMSA primers."
    )


# ============================================================
# IDENTIFY FACING PRIMER PAIR
#
# Expected PCR product is approximately 150 bp.
#
# We allow 80-300 bp initially, then choose:
#
# 1. primers on opposite strands
# 2. facing one another
# 3. best sequence match
# 4. product length closest to 150 bp
# ============================================================

candidates = []


for a in p1:

    for b in p2:

        if a["contig"] != b["contig"]:
            continue

        if a["strand"] == b["strand"]:
            continue


        a_lo = min(a["sstart"], a["send"])
        a_hi = max(a["sstart"], a["send"])

        b_lo = min(b["sstart"], b["send"])
        b_hi = max(b["sstart"], b["send"])


        # ----------------------------------------------------
        # Determine whether primers face each other
        # ----------------------------------------------------

        facing = False
        orientation = None


        if (
            a["strand"] == "plus"
            and b["strand"] == "minus"
            and a_lo < b_hi
        ):

            facing = True
            orientation = "+"


        elif (
            a["strand"] == "minus"
            and b["strand"] == "plus"
            and b_lo < a_hi
        ):

            facing = True
            orientation = "-"


        if not facing:
            continue


        left = min(
            a_lo,
            b_lo
        )

        right = max(
            a_hi,
            b_hi
        )

        product_length = (
            right - left + 1
        )


        if not (
            80 <= product_length <= 300
        ):
            continue


        combined_match = (
            a["length"] +
            b["length"]
        )

        total_mismatch = (
            a["mismatch"] +
            b["mismatch"]
        )

        total_bitscore = (
            a["bitscore"] +
            b["bitscore"]
        )


        # lower tuple = better
        ranking = (
            total_mismatch,
            -combined_match,
            abs(product_length - 150),
            -total_bitscore
        )


        candidates.append(
            (
                ranking,
                a,
                b,
                left,
                right,
                product_length,
                orientation
            )
        )


if not candidates:

    raise RuntimeError(
        "No plausible facing L2256/L2257 primer pair "
        "within 80-300 bp was found."
    )


candidates.sort(
    key=lambda x: x[0]
)


(
    ranking,
    a,
    b,
    left,
    right,
    product_length,
    orientation
) = candidates[0]


contig = a["contig"]


if contig not in genome:

    raise RuntimeError(
        f"Contig {contig} missing from FASTA."
    )


# genomic interval
seq = genome[contig][
    left - 1:right
]


# orient sequence in the L2256 -> L2257 PCR direction
if orientation == "-":
    seq = rc(seq)


# ============================================================
# REPORT
# ============================================================

print()
print("Selected experimental EMSA primer pair")
print("---------------------------------------")

print(
    "L2256:",
    a["contig"],
    a["sstart"],
    a["send"],
    a["strand"]
)

print(
    "L2257:",
    b["contig"],
    b["sstart"],
    b["send"],
    b["strand"]
)

print()
print(
    "Amplicon genomic interval:",
    f"{contig}:{left}-{right}"
)

print(
    "Amplicon orientation:",
    orientation
)

print(
    "Amplicon length:",
    product_length,
    "bp"
)

print()
print("Amplicon sequence:")
print(seq)


# ============================================================
# WRITE FASTA
# ============================================================

with open(
    fasta_out,
    "w"
) as out:

    out.write(
        f">Dsolani_Ds0432_1_Psol_EMSA"
        f"|{contig}:{left}-{right}"
        f"|orientation={orientation}"
        f"|length={product_length}\n"
    )

    for i in range(
        0,
        len(seq),
        70
    ):

        out.write(
            seq[i:i+70] + "\n"
        )


# ============================================================
# WRITE COORDINATES
# ============================================================

with open(
    coord_out,
    "w"
) as out:

    out.write(
        "reference\tcontig\t"
        "amplicon_start\tamplicon_end\t"
        "orientation\tamplicon_length\t"
        "L2256_sstart\tL2256_send\tL2256_strand\t"
        "L2257_sstart\tL2257_send\tL2257_strand\n"
    )

    out.write(
        "Dsolani_Ds0432_1\t"
        f"{contig}\t"
        f"{left}\t"
        f"{right}\t"
        f"{orientation}\t"
        f"{product_length}\t"
        f"{a['sstart']}\t"
        f"{a['send']}\t"
        f"{a['strand']}\t"
        f"{b['sstart']}\t"
        f"{b['send']}\t"
        f"{b['strand']}\n"
    )


print()
print("Written:")
print(fasta_out)
print(coord_out)

PY


    date > "$STATE/04_extract_EMSA_probe.done"

else

    echo
    echo ">>> STEP 04 already completed — skipping"

fi


# ============================================================
# STEP 05
# FUTURE: LOCATE PsolA RELATIVE TO solA
#
# We will add code here AFTER confirming the exact EMSA probe.
# ============================================================

if [[ ! -f "$STATE/05_locate_relative_to_solA.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 05: solA-relative coordinates"
    echo "============================================================"
    echo

    echo "STEP 05 not implemented yet."
    echo "Waiting until exact EMSA fragment is verified."

    # --------------------------------------------------------
    # ADD FUTURE CODE HERE
    # --------------------------------------------------------
    #
    # DO NOT create the checkpoint yet.
    #
    # When this step is implemented successfully, add:
    #
    # date > "$STATE/05_locate_relative_to_solA.done"
    #
    # --------------------------------------------------------

else

    echo
    echo ">>> STEP 05 already completed — skipping"

fi


# ============================================================
# STEP 06
# FUTURE: FIND HOMOLOGOUS Psol REGION IN R. BADENSIS
# ============================================================

if [[ ! -f "$STATE/06_Rouxiella_Psol_homology.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 06: Rouxiella Psol comparison"
    echo "============================================================"
    echo

    echo "STEP 06 not implemented yet."

    # Future:
    #
    # - use exact D. solani EMSA probe
    # - locate homologous region around solA in 20GA0316
    # - repeat across 17 R. badensis genomes
    # - calculate identity / coverage
    #
    # DO NOT touch checkpoint until implemented.

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi


# ============================================================
# CURRENT OUTPUT SUMMARY
# ============================================================

echo
echo "============================================================"
echo " CURRENT ANALYSIS COMPLETE"
echo "============================================================"
echo

echo "Reference genome:"
echo "$GENOME"

echo
echo "Primer hits:"
echo "$PRIMER/Psol_EMSA_primer_hits.tsv"

echo
echo "Exact experimental EMSA probe:"
echo "$PROBE/Psol_EMSA_probe.fna"

echo
echo "Probe coordinates:"
echo "$PROBE/Psol_EMSA_probe_coordinates.tsv"

echo
echo "Finished: $(date)"



# ============================================================
# STEP 05
# SEARCH THE EXACT 150-bp D. SOLANI Psol EMSA FRAGMENT
# ACROSS ALL 17 R. BADENSIS GENOMES
# ============================================================

if [[ ! -f "$STATE/05_Rouxiella_Psol_search.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 05: Searching exact Psol EMSA region in R. badensis"
    echo "============================================================"
    echo


    RB_GENOMES="/scratch/al98750/Roux/02_genomes"
    ACCESSIONS="/scratch/al98750/Roux/01_NCBI/selected_accessions.txt"

    QUERY="$PROBE/Psol_EMSA_probe.fna"

    RAW="$RB/Psol_EMSA_all_hits.tsv"
    BEST="$RESULTS/Psol_EMSA_best_hit_per_genome.tsv"


    # --------------------------------------------------------
    # Header
    # --------------------------------------------------------

    echo -e \
"genome\tstrain_label\tcontig\tpident\talignment_length\tquery_length\tqcov_pct\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstrand\tevalue\tbitscore" \
    > "$BEST"


    while read -r ACC
    do

        [[ -z "$ACC" ]] && continue


        GENOME="$RB_GENOMES/${ACC}.fna"


        if [[ ! -s "$GENOME" ]]
        then

            echo "WARNING: genome not found:"
            echo "$GENOME"

            continue

        fi


        echo
        echo "Searching:"
        echo "$ACC"


        DB="$TMP/${ACC}_Psol_db"

        makeblastdb \
            -in "$GENOME" \
            -dbtype nucl \
            -out "$DB" \
            > /dev/null


        HITFILE="$TMP/${ACC}_Psol_hits.tsv"


        blastn \
            -task blastn \
            -query "$QUERY" \
            -db "$DB" \
            -evalue 10 \
            -word_size 7 \
            -dust no \
            -reward 1 \
            -penalty -2 \
            -gapopen 2 \
            -gapextend 2 \
            -num_threads "$THREADS" \
            -outfmt \
'6 sseqid pident length qlen mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
            > "$HITFILE"


        # ----------------------------------------------------
        # Add accession to all raw results
        # ----------------------------------------------------

        if [[ -s "$HITFILE" ]]
        then

            awk -v G="$ACC" \
                'BEGIN{OFS="\t"} {print G,$0}' \
                "$HITFILE" \
                >> "$RAW"

        fi


        # ----------------------------------------------------
        # Take best hit:
        #
        # prioritize:
        #   1. query coverage
        #   2. bitscore
        #   3. percent identity
        # ----------------------------------------------------

        python - \
            "$ACC" \
            "$HITFILE" \
            "$BEST" <<'PY'

import sys
from pathlib import Path


acc = sys.argv[1]
hitfile = Path(sys.argv[2])
outfile = Path(sys.argv[3])


LABELS = {
    "GCF_020740305.1":
        "R_badensis_20GA0316",

    "GCF_011067265.1":
        "R_badensis_SER3",

    "GCF_002093665.1":
        "R_badensis_DSM100043"
}


label = LABELS.get(
    acc,
    "R_badensis"
)


hits = []


if hitfile.exists():

    for line in hitfile.read_text().splitlines():

        if not line.strip():
            continue

        p = line.split("\t")

        qlen = int(p[3])
        alnlen = int(p[2])

        qcov = (
            100.0 * alnlen / qlen
        )

        hits.append({

            "contig": p[0],

            "pident": float(p[1]),

            "length": alnlen,

            "qlen": qlen,

            "qcov": qcov,

            "mismatch": int(p[4]),

            "gapopen": int(p[5]),

            "qstart": int(p[6]),

            "qend": int(p[7]),

            "sstart": int(p[8]),

            "send": int(p[9]),

            "strand": p[10],

            "evalue": p[11],

            "bitscore": float(p[12])

        })


with open(outfile, "a") as out:

    if not hits:

        out.write(
            f"{acc}\t"
            f"{label}\t"
            f"NO_HIT\t"
            f"NA\tNA\t150\t0\t"
            f"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
        )

    else:

        hits.sort(
            key=lambda h: (
                h["qcov"],
                h["bitscore"],
                h["pident"]
            ),
            reverse=True
        )

        h = hits[0]

        out.write(
            f"{acc}\t"
            f"{label}\t"
            f"{h['contig']}\t"
            f"{h['pident']:.2f}\t"
            f"{h['length']}\t"
            f"{h['qlen']}\t"
            f"{h['qcov']:.2f}\t"
            f"{h['mismatch']}\t"
            f"{h['gapopen']}\t"
            f"{h['qstart']}\t"
            f"{h['qend']}\t"
            f"{h['sstart']}\t"
            f"{h['send']}\t"
            f"{h['strand']}\t"
            f"{h['evalue']}\t"
            f"{h['bitscore']:.1f}\n"
        )

PY


    done < "$ACCESSIONS"


    echo
    echo "============================================================"
    echo "Psol EMSA conservation summary"
    echo "============================================================"
    echo

    column -t -s $'\t' "$BEST"


    date > "$STATE/05_Rouxiella_Psol_search.done"

else

    echo
    echo ">>> STEP 05 already completed — skipping"

fi


# ============================================================
# STEP 06
# LOCATE THE EXPERIMENTAL Psol EMSA PROBE RELATIVE TO solA
# IN D. SOLANI D s0432-1
# ============================================================

if [[ ! -f "$STATE/06_locate_probe_relative_to_solA.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 06: Locating Psol EMSA probe relative to solA"
    echo "============================================================"
    echo


    GB="$REF/Dsolani_Ds0432_1_CP017453.1.gb"

    # --------------------------------------------------------
    # Download annotated GenBank record if needed
    # --------------------------------------------------------

    if [[ ! -s "$GB" ]]
    then

        echo "Downloading annotated CP017453.1 GenBank record..."

        URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP017453.1&rettype=gbwithparts&retmode=text"

        curl \
            -L \
            --fail \
            --retry 3 \
            "$URL" \
            -o "$GB"

    fi


    python - \
        "$GB" \
        "$PROBE/Psol_EMSA_probe_coordinates.tsv" \
        "$RESULTS/Dsolani_probe_nearby_CDS.tsv" \
        "$RESULTS/Dsolani_solA_candidate.tsv" <<'PY'

import sys
import csv
import re
from pathlib import Path


gb_file = Path(sys.argv[1])
probe_file = Path(sys.argv[2])
nearby_out = Path(sys.argv[3])
candidate_out = Path(sys.argv[4])


# ============================================================
# READ EMSA PROBE COORDINATES
# ============================================================

with open(probe_file) as fh:
    probe = next(csv.DictReader(fh, delimiter="\t"))

probe_start = int(probe["amplicon_start"])
probe_end   = int(probe["amplicon_end"])


print()
print("Experimental EMSA probe:")
print(f"{probe_start}-{probe_end}")
print(f"Length: {probe_end-probe_start+1} bp")
print()


# ============================================================
# SIMPLE GENBANK FEATURE PARSER
# ============================================================

features = []

current = None
current_qualifier = None


def finish_feature():

    global current

    if current is not None:
        features.append(current)

    current = None


with open(gb_file) as fh:

    in_features = False

    for raw in fh:

        if raw.startswith("FEATURES"):
            in_features = True
            continue

        if raw.startswith("ORIGIN"):
            finish_feature()
            break

        if not in_features:
            continue


        # ----------------------------------------------------
        # New feature
        # GenBank feature key occupies columns 6-20
        # ----------------------------------------------------

        key = raw[5:21].strip()
        location = raw[21:].strip()


        if key:

            finish_feature()

            current = {
                "type": key,
                "location": location,
                "qualifiers": {}
            }

            current_qualifier = None
            continue


        if current is None:
            continue


        text = raw[21:].strip()


        # ----------------------------------------------------
        # New qualifier
        # ----------------------------------------------------

        if text.startswith("/"):

            if "=" in text:

                q, value = text[1:].split("=", 1)

                value = value.strip().strip('"')

                current["qualifiers"][q] = value

                current_qualifier = q

            else:

                q = text[1:]

                current["qualifiers"][q] = ""

                current_qualifier = q

            continue


        # ----------------------------------------------------
        # Continuation of qualifier
        # ----------------------------------------------------

        if text and current_qualifier:

            old = current["qualifiers"].get(
                current_qualifier,
                ""
            )

            text = text.strip().strip('"')

            current["qualifiers"][current_qualifier] = (
                old + " " + text
            ).strip()


# ============================================================
# PARSE CDS COORDINATES
# ============================================================

cds = []


for f in features:

    if f["type"] != "CDS":
        continue

    loc = f["location"]

    strand = "-1" if loc.startswith("complement") else "+1"


    nums = [
        int(x)
        for x in re.findall(r"\d+", loc)
    ]

    if len(nums) < 2:
        continue


    start = min(nums)
    end = max(nums)

    q = f["qualifiers"]


    cds.append({
        "start": start,
        "end": end,
        "strand": strand,
        "length_nt": end - start + 1,
        "aa_length_est": (end - start + 1) // 3,
        "locus_tag": q.get("locus_tag", ""),
        "gene": q.get("gene", ""),
        "product": q.get("product", "")
    })


# ============================================================
# FIND CDS FEATURES NEAR EMSA PROBE
# ============================================================

WINDOW = 25000

nearby = []


for c in cds:

    if (
        c["end"] >= probe_start - WINDOW
        and
        c["start"] <= probe_end + WINDOW
    ):

        if c["start"] > probe_end:
            distance = c["start"] - probe_end - 1

        elif c["end"] < probe_start:
            distance = probe_start - c["end"] - 1

        else:
            distance = 0


        x = dict(c)

        x["distance_from_probe_bp"] = distance

        nearby.append(x)


nearby.sort(
    key=lambda x: (
        x["start"],
        x["end"]
    )
)


# ============================================================
# WRITE NEARBY CDS TABLE
# ============================================================

fields = [
    "start",
    "end",
    "strand",
    "length_nt",
    "aa_length_est",
    "distance_from_probe_bp",
    "locus_tag",
    "gene",
    "product"
]


with open(
    nearby_out,
    "w",
    newline=""
) as out:

    writer = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(nearby)


print("CDS features within ±25 kb of probe:")
print()

for c in nearby:

    print(
        c["start"],
        c["end"],
        c["strand"],
        f"{c['aa_length_est']} aa",
        f"distance={c['distance_from_probe_bp']} bp",
        c["locus_tag"],
        c["gene"],
        c["product"]
    )


# ============================================================
# IDENTIFY LIKELY solA
#
# SolA is a very large PKS/NRPS protein (~5,000+ aa),
# so the correct CDS should stand out clearly.
# ============================================================

large = [
    c for c in nearby
    if c["aa_length_est"] >= 4000
]


if not large:

    raise RuntimeError(
        "No >4000-aa CDS found near the experimental "
        "Psol EMSA probe. Inspect nearby CDS manually."
    )


large.sort(
    key=lambda x: (
        x["distance_from_probe_bp"],
        -x["aa_length_est"]
    )
)


solA = large[0]


# ============================================================
# CALCULATE PROBE POSITION RELATIVE TO solA START CODON
#
# 0 = first nucleotide of solA CDS
# negative = upstream
# ============================================================

if solA["strand"] == "+1":

    rel1 = probe_start - solA["start"]
    rel2 = probe_end   - solA["start"]

else:

    rel1 = solA["end"] - probe_start
    rel2 = solA["end"] - probe_end


rel_start = min(rel1, rel2)
rel_end   = max(rel1, rel2)


print()
print("============================================================")
print("LIKELY solA")
print("============================================================")
print()

print("Locus tag:", solA["locus_tag"])
print("Gene:", solA["gene"])
print("Product:", solA["product"])

print(
    "CDS:",
    f"{solA['start']}-{solA['end']}"
)

print(
    "Strand:",
    solA["strand"]
)

print(
    "Approximate protein length:",
    solA["aa_length_est"],
    "aa"
)

print(
    "Distance from EMSA probe:",
    solA["distance_from_probe_bp"],
    "bp"
)

print()
print(
    "Experimental probe relative to solA start codon:",
    f"{rel_start} to {rel_end} bp"
)


# ============================================================
# WRITE solA RESULT
# ============================================================

with open(
    candidate_out,
    "w"
) as out:

    out.write(
        "reference\t"
        "locus_tag\t"
        "gene\t"
        "product\t"
        "solA_start\t"
        "solA_end\t"
        "strand\t"
        "protein_length_est_aa\t"
        "probe_start\t"
        "probe_end\t"
        "probe_relative_start_bp\t"
        "probe_relative_end_bp\t"
        "distance_probe_to_CDS_bp\n"
    )

    out.write(
        "Dsolani_Ds0432_1\t"
        f"{solA['locus_tag']}\t"
        f"{solA['gene']}\t"
        f"{solA['product']}\t"
        f"{solA['start']}\t"
        f"{solA['end']}\t"
        f"{solA['strand']}\t"
        f"{solA['aa_length_est']}\t"
        f"{probe_start}\t"
        f"{probe_end}\t"
        f"{rel_start}\t"
        f"{rel_end}\t"
        f"{solA['distance_from_probe_bp']}\n"
    )


print()
print("Written:")
print(nearby_out)
print(candidate_out)

PY


    date > "$STATE/06_locate_probe_relative_to_solA.done"

else

    echo
    echo ">>> STEP 06 already completed — skipping"

fi