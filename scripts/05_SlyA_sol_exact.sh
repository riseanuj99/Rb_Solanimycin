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
# LOCATE EXPERIMENTAL Psol EMSA PROBE RELATIVE TO REAL solA
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
    # Get annotated GenBank record if absent
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
# EMSA PROBE
# ============================================================

with open(probe_file) as fh:
    probe = next(csv.DictReader(fh, delimiter="\t"))

probe_start = int(probe["amplicon_start"])
probe_end   = int(probe["amplicon_end"])


print()
print("Experimental Psol EMSA probe:")
print(f"  {probe_start}-{probe_end}")
print(f"  length = {probe_end-probe_start+1} bp")
print()


# ============================================================
# GENBANK PARSER
#
# IMPORTANT:
# Origin-spanning join()/order() CDS features are skipped.
# The previous parser incorrectly converted such a feature
# into 1..whole_chromosome.
# ============================================================

features = []

current = None
current_qualifier = None
in_features = False


def save_current():

    global current

    if current is not None:
        features.append(current)

    current = None


with open(gb_file) as fh:

    for raw in fh:

        if raw.startswith("FEATURES"):
            in_features = True
            continue

        if raw.startswith("ORIGIN"):
            save_current()
            break

        if not in_features:
            continue


        key = raw[5:21].strip()
        text = raw[21:].rstrip()


        # ----------------------------------------------------
        # New feature
        # ----------------------------------------------------

        if key:

            save_current()

            current = {
                "type": key,
                "location": text.strip(),
                "qualifiers": {}
            }

            current_qualifier = None
            continue


        if current is None:
            continue


        stripped = text.strip()


        # ----------------------------------------------------
        # Qualifier
        # ----------------------------------------------------

        if stripped.startswith("/"):

            if "=" in stripped:

                q, value = stripped[1:].split("=", 1)

                current["qualifiers"][q] = (
                    value.strip().strip('"')
                )

                current_qualifier = q

            else:

                q = stripped[1:]

                current["qualifiers"][q] = ""

                current_qualifier = q

            continue


        # ----------------------------------------------------
        # Location continuation BEFORE qualifiers
        # ----------------------------------------------------

        if current_qualifier is None and stripped:

            current["location"] += stripped
            continue


        # ----------------------------------------------------
        # Qualifier continuation
        # ----------------------------------------------------

        if current_qualifier is not None and stripped:

            old = current["qualifiers"].get(
                current_qualifier,
                ""
            )

            current["qualifiers"][current_qualifier] = (
                old + " " + stripped.strip('"')
            ).strip()


# ============================================================
# PARSE ONLY SIMPLE CONTIGUOUS CDS FEATURES
# ============================================================

cds = []

skipped_complex = 0


for f in features:

    if f["type"] != "CDS":
        continue


    loc = f["location"].replace(" ", "")


    # --------------------------------------------------------
    # Do NOT flatten join/order locations
    # --------------------------------------------------------

    if "join(" in loc or "order(" in loc:

        skipped_complex += 1
        continue


    strand = "-"

    if not loc.startswith("complement("):
        strand = "+"


    # Remove complement wrapper and fuzzy markers
    clean = loc

    if clean.startswith("complement(") and clean.endswith(")"):
        clean = clean[len("complement("):-1]

    clean = clean.replace("<", "").replace(">", "")


    m = re.fullmatch(
        r"(\d+)\.\.(\d+)",
        clean
    )


    if not m:
        continue


    start = int(m.group(1))
    end = int(m.group(2))


    if start > end:
        start, end = end, start


    nt_length = end - start + 1


    # Must resemble a real CDS
    if nt_length <= 0:
        continue


    q = f["qualifiers"]


    cds.append({
        "start": start,
        "end": end,
        "strand": strand,
        "length_nt": nt_length,
        "aa_length_est": nt_length // 3,
        "locus_tag": q.get("locus_tag", ""),
        "gene": q.get("gene", ""),
        "product": q.get("product", "")
    })


print(
    "Complex origin/join CDS features skipped:",
    skipped_complex
)


# ============================================================
# CDS FEATURES NEAR THE EMSA PROBE
# ============================================================

WINDOW = 30000

nearby = []


for c in cds:

    if (
        c["end"] >= probe_start - WINDOW
        and
        c["start"] <= probe_end + WINDOW
    ):

        # genomic distance only
        if c["start"] > probe_end:

            distance = (
                c["start"] -
                probe_end -
                1
            )

        elif c["end"] < probe_start:

            distance = (
                probe_start -
                c["end"] -
                1
            )

        else:

            distance = 0


        # ----------------------------------------------------
        # Is the probe transcriptionally upstream?
        #
        # + gene: promoter must be left of CDS
        # - gene: promoter must be right of CDS
        # ----------------------------------------------------

        if (
            c["strand"] == "+"
            and probe_end < c["start"]
        ):

            upstream_distance = (
                c["start"] -
                probe_end -
                1
            )


        elif (
            c["strand"] == "-"
            and probe_start > c["end"]
        ):

            upstream_distance = (
                probe_start -
                c["end"] -
                1
            )


        else:

            upstream_distance = None


        x = dict(c)

        x["distance_from_probe_bp"] = distance

        x["probe_is_upstream"] = (
            upstream_distance is not None
        )

        x["upstream_distance_bp"] = (
            upstream_distance
            if upstream_distance is not None
            else "NA"
        )


        nearby.append(x)


nearby.sort(
    key=lambda x: (
        x["start"],
        x["end"]
    )
)


# ============================================================
# WRITE NEARBY CDS
# ============================================================

fields = [
    "start",
    "end",
    "strand",
    "length_nt",
    "aa_length_est",
    "distance_from_probe_bp",
    "probe_is_upstream",
    "upstream_distance_bp",
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


print()
print("CDS within +/-30 kb:")
print()


for c in nearby:

    print(
        f"{c['start']}-{c['end']}",
        c["strand"],
        f"{c['aa_length_est']} aa",
        f"upstream={c['probe_is_upstream']}",
        f"upstream_distance={c['upstream_distance_bp']}",
        c["locus_tag"],
        c["gene"],
        c["product"]
    )


# ============================================================
# IDENTIFY solA
#
# Published SolA size ≈ 5358 aa.
#
# Require:
#   4500-6500 aa
#   probe transcriptionally upstream
#   probe <= 5000 bp upstream
# ============================================================

candidates = []


for c in nearby:

    if not (
        4500 <= c["aa_length_est"] <= 6500
    ):
        continue


    if not c["probe_is_upstream"]:
        continue


    d = c["upstream_distance_bp"]


    if d == "NA" or d > 5000:
        continue


    # Prefer:
    # 1. closest to expected 5358 aa
    # 2. closest promoter distance

    score = (
        abs(c["aa_length_est"] - 5358),
        d
    )


    candidates.append(
        (score, c)
    )


if not candidates:

    print()
    print("ERROR:")
    print(
        "No realistic ~5358-aa solA candidate was found "
        "with the EMSA probe upstream."
    )

    print()
    print(
        "This is intentionally stopping rather than "
        "guessing the wrong gene."
    )

    raise SystemExit(1)


candidates.sort(
    key=lambda x: x[0]
)


solA = candidates[0][1]


# ============================================================
# POSITION OF EMSA PROBE RELATIVE TO solA START CODON
#
# CDS first nucleotide = position 0
# upstream positions = negative
# ============================================================

if solA["strand"] == "+":

    rel1 = probe_start - solA["start"]
    rel2 = probe_end   - solA["start"]


else:

    # On minus strand, start codon is at genomic CDS end.
    rel1 = solA["end"] - probe_start
    rel2 = solA["end"] - probe_end


rel_start = min(rel1, rel2)
rel_end   = max(rel1, rel2)


# ============================================================
# SANITY CHECK
# ============================================================

if rel_end >= 0:

    raise RuntimeError(
        "Probe is not fully upstream of selected solA."
    )


print()
print("============================================================")
print("REALISTIC solA CANDIDATE")
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
    "Protein length:",
    solA["aa_length_est"],
    "aa"
)

print(
    "Probe-to-CDS distance:",
    solA["upstream_distance_bp"],
    "bp"
)

print()
print(
    "150-bp EMSA probe relative to solA ATG:",
    f"{rel_start} to {rel_end}"
)


# ============================================================
# WRITE RESULT
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
        f"{solA['upstream_distance_bp']}\n"
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

# ============================================================


# ============================================================
# STEP 07
# EXTRACT SolA UPSTREAM REGIONS USING VERIFIED sol LOCI
# FROM THE EXISTING sol-HGT ANALYSIS
# ============================================================

if [[ ! -f "$STATE/07_Rouxiella_solA_upstream.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 07: Extracting upstream DNA from verified solA loci"
    echo "============================================================"
    echo


    SOL_HITS="/scratch/al98750/Roux/07_sol_HGT/05_results/sol_gene_hits.tsv"
    SOL_LOCI="/scratch/al98750/Roux/07_sol_HGT/05_results/sol_loci.tsv"

    RB_GENOMES="/scratch/al98750/Roux/02_genomes"
    ACCESSIONS="/scratch/al98750/Roux/01_NCBI/selected_accessions.txt"

    UPSTREAM_FASTA="$RB/Rouxiella_SolA_upstream600.fna"

    UPSTREAM_META="$RESULTS/Rouxiella_SolA_upstream600_metadata.tsv"

    UPSTREAM=600


    python - \
        "$SOL_HITS" \
        "$SOL_LOCI" \
        "$RB_GENOMES" \
        "$ACCESSIONS" \
        "$UPSTREAM_FASTA" \
        "$UPSTREAM_META" \
        "$UPSTREAM" <<'PY'

import sys
import csv
from pathlib import Path


sol_hits_file = Path(sys.argv[1])
sol_loci_file = Path(sys.argv[2])
genome_root = Path(sys.argv[3])
accessions_file = Path(sys.argv[4])

fasta_out = Path(sys.argv[5])
meta_out = Path(sys.argv[6])

UPSTREAM = int(sys.argv[7])


# ============================================================
# STRAIN LABELS
# ============================================================

LABELS = {

    "GCF_020740205.1":
        "R_badensis_20GA0206",

    "GCF_020740305.1":
        "R_badensis_20GA0316",

    "GCF_020740185.1":
        "R_badensis_20GA0366",

    "GCF_020739655.1":
        "R_badensis_20GA0383",

    "GCF_020740225.1":
        "R_badensis_20GA0386",

    "GCF_020739645.1":
        "R_badensis_21GA0654",

    "GCF_020740165.1":
        "R_badensis_Rb14-3",

    "GCF_002093665.1":
        "R_badensis_DSM100043",

    "GCF_011067265.1":
        "R_badensis_SER3"
}


# ============================================================
# HELPERS
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
            "ACGTNacgtn",
            "TGCANtgcan"
        )
    )[::-1]


# ============================================================
# SELECTED GENOMES
# ============================================================

selected = {

    x.strip()

    for x in accessions_file.read_text().splitlines()

    if x.strip()

}


# ============================================================
# READ sol CLUSTER STATUS
# ============================================================

cluster_status = {}


with open(sol_loci_file) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )


    for r in reader:

        cluster_status[r["accession"]] = {

            "status":
                r["status"],

            "genes_detected":
                r["genes_detected"],

            "cluster_contig":
                r["contig"],

            "cluster_start_1based":
                r["start_1based"],

            "cluster_end":
                r["end"]

        }


# ============================================================
# READ VERIFIED solA HITS
#
# IMPORTANT:
# sol_gene_hits coordinates are 0-based start / end-exclusive
# style coordinates.
#
# Evidence:
# sol_loci reports start_1based = solA start + 1.
#
# Thus:
#
#   start0 = start
#   CDS in Python = sequence[start:end]
#
# ============================================================

solA_hits = {}


with open(sol_hits_file) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )


    for r in reader:

        if r["gene"] != "solA":
            continue

        acc = r["accession"]


        if acc not in selected:
            continue


        solA_hits[acc] = r


# ============================================================
# EXTRACT UPSTREAM REGION
# ============================================================

rows = []
sequences = []


for acc in sorted(selected):

    if acc not in solA_hits:

        print(
            "No verified solA hit:",
            acc
        )

        continue


    h = solA_hits[acc]


    contig = h["contig"]

    start0 = int(h["start"])

    end0 = int(h["end"])

    strand_raw = h["strand"]


    if strand_raw in {"1", "+", "+1"}:
        strand = "+"

    elif strand_raw in {"-1", "-"}:
        strand = "-"

    else:

        print(
            "WARNING: unexpected strand:",
            acc,
            strand_raw
        )

        continue


    # --------------------------------------------------------
    # Locate genome FASTA
    # --------------------------------------------------------

    genome_file = (
        genome_root /
        f"{acc}.fna"
    )


    if not genome_file.exists():

        matches = list(
            genome_root.glob(
                f"{acc}*.fna"
            )
        )


        if not matches:

            print(
                "WARNING: genome FASTA missing:",
                acc
            )

            continue


        genome_file = matches[0]


    genome = read_fasta(
        genome_file
    )


    if contig not in genome:

        print(
            "WARNING: contig missing:",
            acc,
            contig
        )

        continue


    chrom = genome[contig]


    # --------------------------------------------------------
    # Extract exactly 600 nt upstream of the transcriptional
    # start side of SolA.
    #
    # sol_gene_hits:
    #
    #   start0 = 0-based CDS start
    #   end0   = end-exclusive CDS end
    #
    # + strand:
    #    promoter = start0-600 : start0
    #
    # - strand:
    #    promoter = end0 : end0+600
    #    then reverse complement
    #
    # Output orientation is always:
    #
    #    -600 ............. -1 | solA
    #
    # --------------------------------------------------------

    if strand == "+":

        left0 = max(
            0,
            start0 - UPSTREAM
        )

        right0 = start0

        seq = chrom[
            left0:right0
        ]


        genomic_left_1based = (
            left0 + 1
        )

        genomic_right_1based = (
            right0
        )


        solA_start_1based = (
            start0 + 1
        )

        solA_end_1based = (
            end0
        )


    else:

        left0 = end0

        right0 = min(
            len(chrom),
            end0 + UPSTREAM
        )


        seq = rc(
            chrom[
                left0:right0
            ]
        )


        genomic_left_1based = (
            left0 + 1
        )

        genomic_right_1based = (
            right0
        )


        solA_start_1based = (
            start0 + 1
        )

        solA_end_1based = (
            end0
        )


    # --------------------------------------------------------
    # Cluster information
    # --------------------------------------------------------

    cluster = cluster_status.get(
        acc,
        {}
    )


    label = LABELS.get(
        acc,
        "R_badensis"
    )


    sequences.append(
        (
            acc,
            label,
            seq
        )
    )


    rows.append({

        "genome":
            acc,

        "strain_label":
            label,

        "cluster_status":
            cluster.get(
                "status",
                "NA"
            ),

        "sol_genes_detected":
            cluster.get(
                "genes_detected",
                "NA"
            ),

        "solA_protein":
            h["protein"],

        "solA_contig":
            contig,

        "solA_start_1based":
            solA_start_1based,

        "solA_end_1based":
            solA_end_1based,

        "solA_strand":
            strand,

        "solA_identity_pct":
            h["identity"],

        "solA_coverage_pct":
            h["coverage"],

        "solA_evalue":
            h["evalue"],

        "upstream_length":
            len(seq),

        "upstream_genomic_left_1based":
            genomic_left_1based,

        "upstream_genomic_right_1based":
            genomic_right_1based

    })


# ============================================================
# WRITE FASTA
# ============================================================

with open(
    fasta_out,
    "w"
) as out:

    for acc, label, seq in sequences:

        out.write(
            f">{acc}|{label}\n"
        )


        for i in range(
            0,
            len(seq),
            70
        ):

            out.write(
                seq[i:i+70]
                + "\n"
            )


# ============================================================
# WRITE METADATA
# ============================================================

fields = [

    "genome",
    "strain_label",
    "cluster_status",
    "sol_genes_detected",
    "solA_protein",
    "solA_contig",
    "solA_start_1based",
    "solA_end_1based",
    "solA_strand",
    "solA_identity_pct",
    "solA_coverage_pct",
    "solA_evalue",
    "upstream_length",
    "upstream_genomic_left_1based",
    "upstream_genomic_right_1based"

]


with open(
    meta_out,
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
print(
    "Verified SolA upstream regions extracted:",
    len(rows)
)

print()
print("Written:")
print(fasta_out)
print(meta_out)

PY


    echo
    echo "============================================================"
    echo "Verified Rouxiella SolA upstream regions"
    echo "============================================================"
    echo

    column -t -s $'\t' \
        "$UPSTREAM_META"


    date > "$STATE/07_Rouxiella_solA_upstream.done"

else

    echo
    echo ">>> STEP 07 already completed — skipping"

fi

# ============================================================
# STEP 08
# COMPARE THE POSITION-EQUIVALENT SlyA-BOUND PsolA REGION
#
# Experimental D. solani EMSA probe:
#   relative to SolA ATG = -436 to -287
#   length = 150 bp
#
# Question:
# What occupies the same SolA-relative position in R. badensis?
# ============================================================

if [[ ! -f "$STATE/08_Psol_position_equivalent.done" ]]
then

    echo
    echo "============================================================"
    echo "STEP 08: Comparing homologous-position PsolA regions"
    echo "============================================================"
    echo


    RB_UPSTREAM="$RB/Rouxiella_SolA_upstream600.fna"

    DSOLA_TABLE="$RESULTS/Dsolani_solA_candidate.tsv"

    DS_UPSTREAM="$RB/Dsolani_SolA_upstream600.fna"

    COMBINED_600="$RB/SolA_upstream600_all.fna"

    SITE_FASTA="$RB/Psol_position_equivalent_150bp.fna"

    SITE_ALN="$RB/Psol_position_equivalent_150bp.mafft.fna"

    SITE_RESULTS="$RESULTS/Psol_position_equivalent_conservation.tsv"


    # ========================================================
    # 08A
    # EXTRACT D. SOLANI 600-bp UPSTREAM REGION
    # AND THE EXACT -436 TO -287 INTERVAL
    # ========================================================

    python - \
        "$GENOME" \
        "$DSOLA_TABLE" \
        "$PROBE/Psol_EMSA_probe.fna" \
        "$RB_UPSTREAM" \
        "$DS_UPSTREAM" \
        "$COMBINED_600" \
        "$SITE_FASTA" <<'PY'

import sys
import csv
from pathlib import Path


genome_file = Path(sys.argv[1])
sola_table = Path(sys.argv[2])
probe_fasta = Path(sys.argv[3])
rb_upstream_file = Path(sys.argv[4])

ds_upstream_out = Path(sys.argv[5])
combined_out = Path(sys.argv[6])
site_out = Path(sys.argv[7])


UPSTREAM = 600


# ============================================================
# HELPERS
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
            "ACGTNacgtn",
            "TGCANtgcan"
        )
    )[::-1]


# ============================================================
# READ D. SOLANI SolA INFORMATION
# ============================================================

with open(sola_table) as fh:

    solA = next(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


solA_start = int(
    solA["solA_start"]
)

solA_end = int(
    solA["solA_end"]
)

strand = solA["strand"]


probe_rel_start = int(
    solA["probe_relative_start_bp"]
)

probe_rel_end = int(
    solA["probe_relative_end_bp"]
)


print()
print("D. solani experimental Psol region:")
print(
    f"{probe_rel_start} to "
    f"{probe_rel_end} relative to SolA ATG"
)


# ============================================================
# READ D. SOLANI CHROMOSOME
# ============================================================

genome = read_fasta(
    genome_file
)


if len(genome) != 1:

    raise RuntimeError(
        "Expected one D. solani chromosome."
    )


chrom = next(
    iter(genome.values())
)


# ============================================================
# EXTRACT 600 bp UPSTREAM
# ============================================================

if strand in {"+", "+1"}:

    start0 = solA_start - 1

    left0 = start0 - UPSTREAM

    if left0 < 0:
        raise RuntimeError(
            "D. solani SolA too close to chromosome start."
        )

    ds_upstream = chrom[
        left0:start0
    ]


elif strand in {"-", "-1"}:

    # First nucleotide of minus-strand CDS is genomic solA_end

    left0 = solA_end

    right0 = solA_end + UPSTREAM

    ds_upstream = rc(
        chrom[
            left0:right0
        ]
    )


else:

    raise RuntimeError(
        f"Unexpected SolA strand: {strand}"
    )


if len(ds_upstream) != 600:

    raise RuntimeError(
        f"D. solani upstream length = {len(ds_upstream)}, "
        "expected 600."
    )


# ============================================================
# CONVERT RELATIVE COORDINATES INTO POSITIONS
# WITHIN THE 600-bp UPSTREAM SEQUENCE
#
# Sequence orientation:
#
# position 1                 position 600
#   -600  .......................  -1
#                                  | ATG
#
# For example:
#
# -436 -> sequence position 165
# -287 -> sequence position 314
# ============================================================

site_start_1based = (
    UPSTREAM
    + probe_rel_start
    + 1
)

site_end_1based = (
    UPSTREAM
    + probe_rel_end
    + 1
)


site_start0 = (
    site_start_1based - 1
)

site_end0 = (
    site_end_1based
)


print()
print(
    "Equivalent positions inside 600-bp upstream sequence:",
    f"{site_start_1based}-{site_end_1based}"
)


ds_site = ds_upstream[
    site_start0:site_end0
]


print(
    "Extracted D. solani site length:",
    len(ds_site)
)


if len(ds_site) != 150:

    raise RuntimeError(
        f"Expected 150 bp, got {len(ds_site)}"
    )


# ============================================================
# VERIFY AGAINST THE ACTUAL EMSA PCR AMPLICON
# ============================================================

probe_records = read_fasta(
    probe_fasta
)


experimental_probe = next(
    iter(probe_records.values())
)


if ds_site != experimental_probe:

    print()
    print("ERROR:")
    print(
        "The -436/-287 sequence does not exactly match "
        "the reconstructed EMSA probe."
    )

    print()
    print("Extracted:")
    print(ds_site)

    print()
    print("Experimental:")
    print(experimental_probe)

    raise SystemExit(1)


print()
print(
    "PASS: extracted -436/-287 region exactly matches "
    "the 150-bp experimental EMSA probe."
)


# ============================================================
# WRITE D. SOLANI 600 bp
# ============================================================

with open(
    ds_upstream_out,
    "w"
) as out:

    out.write(
        ">Dsolani_Ds0432_1|experimental_reference\n"
    )

    for i in range(
        0,
        len(ds_upstream),
        70
    ):

        out.write(
            ds_upstream[i:i+70]
            + "\n"
        )


# ============================================================
# READ ROUXIELLA 600-bp UPSTREAM SEQUENCES
# ============================================================

rb = read_fasta(
    rb_upstream_file
)


# ============================================================
# WRITE COMBINED 600-bp FASTA
# ============================================================

with open(
    combined_out,
    "w"
) as out:

    out.write(
        ">Dsolani_Ds0432_1|experimental_reference\n"
        + ds_upstream
        + "\n"
    )

    for header, seq in rb.items():

        if len(seq) != 600:

            raise RuntimeError(
                f"{header}: upstream length "
                f"{len(seq)}, expected 600."
            )

        out.write(
            f">{header}\n{seq}\n"
        )


# ============================================================
# EXTRACT EXACT POSITION-EQUIVALENT 150-bp REGION
# FROM EVERY ROUXIELLA UPSTREAM SEQUENCE
#
# This is positional homology:
#
#       -436 ............... -287
#
# relative to each strain's own SolA ATG.
# ============================================================

with open(
    site_out,
    "w"
) as out:

    out.write(
        ">Dsolani_Ds0432_1|experimental_reference\n"
        + ds_site
        + "\n"
    )


    for header, seq in rb.items():

        site = seq[
            site_start0:site_end0
        ]


        if len(site) != 150:

            raise RuntimeError(
                f"{header}: site length "
                f"{len(site)}, expected 150."
            )


        out.write(
            f">{header}\n"
            f"{site}\n"
        )


print()
print("Written:")
print(ds_upstream_out)
print(combined_out)
print(site_out)

PY


    # ========================================================
    # 08B
    # ALIGN THE POSITION-EQUIVALENT 150-bp REGIONS
    # ========================================================

    echo
    echo "Aligning 150-bp position-equivalent regions..."


    mafft \
        --localpair \
        --maxiterate 1000 \
        --thread "$THREADS" \
        "$SITE_FASTA" \
        > "$SITE_ALN"


    # ========================================================
    # 08C
    # CALCULATE CONSERVATION TO D. SOLANI
    # AND TO THE FOCAL ONION STRAIN 20GA0316
    # ========================================================

    python - \
        "$SITE_ALN" \
        "$SITE_RESULTS" <<'PY'

import sys
import csv
from pathlib import Path


aln_file = Path(sys.argv[1])
out_file = Path(sys.argv[2])


DS_REF = "Dsolani_Ds0432_1"

RB_REF = "GCF_020740305.1"


# ============================================================
# BIOLOGICAL GROUPS
# ============================================================

ONION_GA = {

    "GCF_020740205.1": "20GA0206",
    "GCF_020740305.1": "20GA0316",
    "GCF_020740185.1": "20GA0366",
    "GCF_020739655.1": "20GA0383",
    "GCF_020740225.1": "20GA0386",
    "GCF_020739645.1": "21GA0654"

}


COMPARATORS = {

    "GCF_020740165.1":
        "Rb14-3",

    "GCF_002093665.1":
        "DSM100043",

    "GCF_011067265.1":
        "SER3"

}


# ============================================================
# FASTA READER
# ============================================================

def read_fasta(path):

    d = {}

    name = None
    seq = []


    for line in path.read_text().splitlines():

        line = line.strip()

        if not line:
            continue


        if line.startswith(">"):

            if name is not None:
                d[name] = "".join(seq).upper()

            name = line[1:].split()[0]
            seq = []

        else:

            seq.append(line)


    if name is not None:
        d[name] = "".join(seq).upper()


    return d


aln = read_fasta(
    aln_file
)


# ============================================================
# IDENTIFY ACCESSION FROM HEADER
# ============================================================

def accession(header):

    return header.split("|")[0]


ds_header = next(
    h for h in aln
    if accession(h) == DS_REF
)


rb_header = next(
    h for h in aln
    if accession(h) == RB_REF
)


ds = aln[ds_header]
rb = aln[rb_header]


# ============================================================
# PAIRWISE CONSERVATION FUNCTION
# ============================================================

def compare(reference, query):

    cols = [
        i
        for i, nt in enumerate(reference)
        if nt != "-"
    ]


    matches = 0
    nongap = 0


    for i in cols:

        q = query[i]
        r = reference[i]


        if q != "-":
            nongap += 1


        if q != "-" and q == r:
            matches += 1


    identity = (
        100.0 * matches / len(cols)
    )

    coverage = (
        100.0 * nongap / len(cols)
    )


    return (
        identity,
        coverage
    )


# ============================================================
# BUILD TABLE
# ============================================================

rows = []


for header, seq in aln.items():

    acc = accession(header)


    if acc == DS_REF:

        strain = "D_s0432-1"
        group = "D_solani_reference"

    elif acc in ONION_GA:

        strain = ONION_GA[acc]
        group = "GA_onion"

    elif acc in COMPARATORS:

        strain = COMPARATORS[acc]
        group = "R_badensis_comparator"

    else:

        strain = acc
        group = "other_R_badensis"


    ds_id, ds_cov = compare(
        ds,
        seq
    )


    rb_id, rb_cov = compare(
        rb,
        seq
    )


    ungapped = seq.replace(
        "-",
        ""
    )


    rows.append({

        "genome":
            acc,

        "strain":
            strain,

        "group":
            group,

        "identity_to_Dsolani_EMSA_region_pct":
            round(ds_id, 2),

        "coverage_vs_Dsolani_EMSA_region_pct":
            round(ds_cov, 2),

        "identity_to_20GA0316_region_pct":
            round(rb_id, 2),

        "coverage_vs_20GA0316_region_pct":
            round(rb_cov, 2),

        "ungapped_length":
            len(ungapped),

        "position_equivalent_sequence":
            ungapped

    })


# Sort biologically
group_order = {

    "D_solani_reference": 0,
    "GA_onion": 1,
    "R_badensis_comparator": 2,
    "other_R_badensis": 3

}


rows.sort(
    key=lambda r: (
        group_order.get(
            r["group"],
            99
        ),
        r["strain"]
    )
)


fields = [

    "genome",
    "strain",
    "group",
    "identity_to_Dsolani_EMSA_region_pct",
    "coverage_vs_Dsolani_EMSA_region_pct",
    "identity_to_20GA0316_region_pct",
    "coverage_vs_20GA0316_region_pct",
    "ungapped_length",
    "position_equivalent_sequence"

]


with open(
    out_file,
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
print("Written:")
print(out_file)

PY


    echo
    echo "============================================================"
    echo "Psol position-equivalent conservation"
    echo "============================================================"
    echo

    column -t -s $'\t' \
        "$SITE_RESULTS"


    date > "$STATE/08_Psol_position_equivalent.done"

else

    echo
    echo ">>> STEP 08 already completed — skipping"

fi
