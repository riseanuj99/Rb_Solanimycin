#!/usr/bin/env bash
# Expanded sol-cluster taxon discovery (NCBI remote BLASTP)
# Usage:
#   bash 11_expand_sol_taxa.sh pilot     # solC only
#   bash 11_expand_sol_taxa.sh all       # six markers, resumes pilot
# Optional: SEARCH_SCOPE=allbacteria bash 11_expand_sol_taxa.sh all
# This is DISCOVERY ONLY: multi-marker hits do not establish a linked sol cluster.

set -Eeuo pipefail

MODE="${1:-all}"
case "$MODE" in
  pilot) GENES=(solC) ;;
  all)   GENES=(solC solE solL solF solH solA) ;;
  *) echo "Usage: bash $0 [pilot|all]" >&2; exit 2 ;;
esac

HGT="${HGT:-/scratch/al98750/Roux/07_sol_HGT}"
REF="$HGT/11_expanded_taxa/00_reference/sol_all13_reference.faa"
SCOPE="${SEARCH_SCOPE:-nonrouxiella}"
case "$SCOPE" in
  nonrouxiella) ENTREZ_QUERY='Bacteria[Organism] NOT Rouxiella[Organism]' ;;
  allbacteria)  ENTREZ_QUERY='Bacteria[Organism]' ;;
  *) echo "SEARCH_SCOPE must be nonrouxiella or allbacteria" >&2; exit 2 ;;
esac

OUT="$HGT/11_expanded_taxa/01_discovery/$SCOPE"
QUERY_DIR="$HGT/11_expanded_taxa/00_reference/individual_queries"
MAX_TARGETS="${MAX_TARGETS:-500}"
TIME_LIMIT="${TIME_LIMIT:-5400}"          # seconds per attempt (90 minutes)
MAX_ATTEMPTS="${MAX_ATTEMPTS:-2}"
WAIT_SECONDS="${WAIT_SECONDS:-40}"       # polite gap between NCBI requests
OUTFMT='6 qseqid saccver pident length qlen slen qcovhsp evalue bitscore staxids sscinames stitle'

for tool in blastp python timeout; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Missing required tool: $tool" >&2; exit 1; }
done
[[ -s "$REF" ]] || { echo "Missing reference: $REF" >&2; exit 1; }
mkdir -p "$OUT" "$QUERY_DIR"

# Extract six exact full-length reference proteins (not trimmed alignments).
python - "$REF" "$QUERY_DIR" <<'PY'
from pathlib import Path
import sys
from Bio import SeqIO

source, target = Path(sys.argv[1]), Path(sys.argv[2])
ids = ('solC', 'solE', 'solL', 'solF', 'solH', 'solA')
records = {}
for record in SeqIO.parse(source, 'fasta'):
    if record.id in records:
        raise SystemExit(f'Duplicate reference protein ID: {record.id}')
    records[record.id] = record
missing = [g for g in ids if g not in records]
if missing:
    raise SystemExit('Missing reference proteins: ' + ', '.join(missing))
for gene in ids:
    path = target / f'{gene}.faa'
    SeqIO.write([records[gene]], path, 'fasta')
    print(f'QUERY\t{gene}\t{len(records[gene].seq)} aa\t{path}')
PY

printf '\nSEARCH MODE: %s\nSCOPE: %s\nNCBI database: nr\nOUTPUT: %s\n' "$MODE" "$SCOPE" "$OUT"

for gene in "${GENES[@]}"; do
  output="$OUT/${gene}_nr.tsv"
  donefile="$OUT/${gene}.done"
  logfile="$OUT/${gene}.log"
  if [[ -f "$donefile" && -f "$output" ]]; then
    printf '\nSKIP %s: completed (%s hits)\n' "$gene" "$(wc -l < "$output")"
    continue
  fi

  success=0
  for ((attempt=1; attempt<=MAX_ATTEMPTS; attempt++)); do
    printf '\nSTART %s (attempt %d/%d): %s\n' "$gene" "$attempt" "$MAX_ATTEMPTS" "$(date)"
    tmp="$OUT/${gene}.partial.tsv"
    rm -f "$tmp"
    if timeout "${TIME_LIMIT}s" blastp \
        -remote -db nr \
        -query "$QUERY_DIR/${gene}.faa" \
        -entrez_query "$ENTREZ_QUERY" \
        -evalue 1e-10 \
        -max_target_seqs "$MAX_TARGETS" \
        -outfmt "$OUTFMT" \
        -out "$tmp" >"$logfile" 2>&1; then
      mv -f "$tmp" "$output"
      touch "$donefile"
      printf 'PASS %s: %s hits (%s)\n' "$gene" "$(wc -l < "$output")" "$(date)"
      success=1
      break
    else
      rc=$?
      rm -f "$tmp"
      printf 'FAILED %s attempt %d (exit %s); log: %s\n' "$gene" "$attempt" "$rc" "$logfile" >&2
      tail -n 8 "$logfile" >&2 || true
      if ((attempt < MAX_ATTEMPTS)); then sleep "$WAIT_SECONDS"; fi
    fi
  done
  if ((success == 0)); then
    echo "WARNING: $gene incomplete. Re-running this script will retry it." >&2
  fi
  sleep "$WAIT_SECONDS"
done

# Aggregate completed searches. Species-level overlap is only a lead, NOT proof
# that different protein hits belong to the same genome or locus.
python - "$OUT" <<'PY'
from collections import defaultdict
from pathlib import Path
import csv, sys

out = Path(sys.argv[1])
columns = ('gene','qseqid','saccver','pident','length','qlen','slen',
           'qcovhsp','evalue','bitscore','staxids','sscinames','stitle')
combined = out / 'all_discovery_hits.tsv'
summary = out / 'species_marker_screen_NOT_cluster_calls.tsv'
completed = []
by_species = defaultdict(lambda: {'genes':set(), 'proteins':set()})

with combined.open('w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(columns)
    for gene in ('solC','solE','solL','solF','solH','solA'):
        path = out / f'{gene}_nr.tsv'
        if not path.is_file() or not (out / f'{gene}.done').is_file():
            continue
        completed.append(gene)
        with path.open(newline='') as inp:
            for row in csv.reader(inp, delimiter='\t'):
                if len(row) != 12:
                    continue
                writer.writerow([gene, *row])
                species = row[10].strip() or 'UNKNOWN'
                by_species[species]['genes'].add(gene)
                by_species[species]['proteins'].add(row[1])

with summary.open('w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(('reported_taxon','number_of_markers','markers',
                     'distinct_protein_accessions','interpretation'))
    for species, v in sorted(by_species.items(),
                             key=lambda x: (-len(x[1]['genes']),
                                            -len(x[1]['proteins']), x[0])):
        writer.writerow((species, len(v['genes']), ','.join(sorted(v['genes'])),
                         len(v['proteins']), 'taxon-level screen only; NOT cluster verified'))

print('\n===== SEARCH SUMMARY =====')
print('Completed markers:', ', '.join(completed) or 'none')
print('Reported taxa:', len(by_species))
print('Combined hits:', combined)
print('Taxon-level screening:', summary)
print('WARNING: Multi-marker taxon hits do NOT prove same genome or linked sol locus.')
PY

echo
echo '===== RESULT COUNTS ====='
for gene in "${GENES[@]}"; do
  if [[ -f "$OUT/${gene}.done" ]]; then
    printf '%-6s %s hits\n' "$gene" "$(wc -l < "$OUT/${gene}_nr.tsv")"
  else
    printf '%-6s INCOMPLETE (see %s/%s.log)\n' "$gene" "$OUT" "$gene"
  fi
done

echo "FINISHED: $(date)"
