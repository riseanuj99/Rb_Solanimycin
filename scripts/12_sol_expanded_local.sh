#!/usr/bin/env bash
# Expanded sol-cluster discovery, NCBI Datasets -> local DIAMOND -> genomic-linkage screen.
# Modes: plan | download | analyze
# IMPORTANT: Run analyze on a Sapelo2 COMPUTE node, not a login node.
set -euo pipefail

MODE="${1:-}"
if [[ ! "$MODE" =~ ^(plan|download|analyze)$ ]]; then
  echo 'Usage: bash 12_sol_expanded_local.sh {plan|download|analyze}' >&2
  exit 2
fi

ROOT="${SOL_ROOT:-/scratch/al98750/Roux/07_sol_HGT/11_expanded_taxa/02_local_discovery}"
REF="${SOL_REF:-/scratch/al98750/Roux/07_sol_HGT/11_expanded_taxa/00_reference/sol_all13_reference.faa}"
OLD_MANIFEST="${SOL_OLD_MANIFEST:-/scratch/al98750/Roux/07_sol_HGT/05_results/genome_manifest.tsv}"
PER_GENUS="${SOL_PER_GENUS:-12}"
PER_SPECIES="${SOL_PER_SPECIES:-2}"
META_LIMIT="${SOL_META_LIMIT:-300}"
BATCH_SIZE="${SOL_BATCH_SIZE:-20}"
THREADS="${SLURM_CPUS_PER_TASK:-${SOL_THREADS:-8}}"
# Taxonomically broad *pilot*, NOT a comprehensive bacterial genome search.
GENERA="${SOL_GENERA:-Rouxiella Pantoea Erwinia Dickeya Pectobacterium Serratia Yersinia Rahnella Hafnia Tatumella Mixta Lonsdalea Brenneria Rosenbergiella Cedecea Enterobacter Klebsiella Citrobacter Kosakonia Kluyvera}"

META="$ROOT/00_metadata"
BATCH="$ROOT/01_batches"
DATA="$ROOT/02_downloads"
ANALYSIS="$ROOT/03_analysis"
mkdir -p "$META" "$BATCH" "$DATA" "$ANALYSIS"

require() { command -v "$1" >/dev/null || { echo "Missing command: $1" >&2; exit 1; }; }

if [[ "$MODE" == plan ]]; then
  require datasets; require python3; require timeout
  echo '== Collecting NCBI RefSeq metadata =='
  for GENUS in $GENERA; do
    FILE="$META/${GENUS}.jsonl"
    if [[ -s "$FILE" ]]; then echo "CACHED $GENUS"; continue; fi
    echo "FETCH $GENUS (limit=$META_LIMIT)"
    for TRY in 1 2; do
      if timeout 150s datasets summary genome taxon "$GENUS" \
          --assembly-source RefSeq --annotated --limit "$META_LIMIT" \
          --as-json-lines > "${FILE}.tmp" 2> "$META/${GENUS}.err"; then
        if [[ -s "${FILE}.tmp" ]]; then mv "${FILE}.tmp" "$FILE"; break; fi
      fi
      rm -f "${FILE}.tmp"
      echo "  attempt $TRY failed: $GENUS" >&2
      sleep 3
    done
    [[ -s "$FILE" ]] || echo "WARNING: no metadata for $GENUS" >&2
    sleep 1
  done
  python3 - "$META" "$OLD_MANIFEST" "$PER_GENUS" "$PER_SPECIES" $GENERA <<'PY'
import csv, hashlib, json, pathlib, re, sys
meta=pathlib.Path(sys.argv[1]); old=pathlib.Path(sys.argv[2]); cap=int(sys.argv[3]); per_sp=int(sys.argv[4]); genera=sys.argv[5:]
def hashkey(s): return hashlib.sha256(s.encode()).hexdigest()
def quality(d):
    level=(d.get('assembly_info') or {}).get('assembly_level','')
    return ({'Complete Genome':0,'Chromosome':1,'Scaffold':2,'Contig':3}.get(level,4),hashkey(d.get('accession','')))
selected={}; stats=[]
for genus in genera:
    path=meta/f'{genus}.jsonl'
    if not path.exists(): stats.append((genus,0,0)); continue
    rows=[]
    for line in path.open():
        try: d=json.loads(line)
        except json.JSONDecodeError: continue
        acc=d.get('accession',''); org=(d.get('organism') or {}).get('organism_name','')
        if not acc.startswith(('GCF_','GCA_')) or not (org==genus or org.startswith(genus+' ')): continue
        rows.append(d)
    species={}
    for d in rows:
        species.setdefault((d.get('organism') or {}).get('organism_name','unknown'),[]).append(d)
    for key in species: species[key].sort(key=quality)
    # Round-robin species so a genus is not dominated by one sequenced species.
    sp_order=sorted(species,key=hashkey); picked=[]
    for i in range(per_sp):
        for sp in sp_order:
            if len(picked)>=cap: break
            if len(species[sp])>i: picked.append(species[sp][i])
        if len(picked)>=cap: break
    for d in picked:
        acc=d['accession']; info=d.get('assembly_info') or {}
        selected.setdefault(acc,(genus,(d.get('organism') or {}).get('organism_name',''),info.get('assembly_level',''),'taxon_sample'))
    stats.append((genus,len(rows),len(picked)))
# Retain existing HGT-panel accessions as controls, even if not in the new sample.
if old.exists():
    with old.open(newline='') as f:
        rd=csv.DictReader(f,delimiter='\t')
        if rd.fieldnames:
            cols=[c for c in rd.fieldnames if c and ('accession' in c.lower() or c.lower()=='genome')]
            for row in rd:
                for col in cols:
                    acc=(row.get(col) or '').strip()
                    if re.fullmatch(r'GC[AF]_\d+\.\d+',acc):
                        selected.setdefault(acc,('existing_panel','','','existing_HGT_panel'))
                        break
with (meta/'selected_genomes.tsv').open('w') as f:
    f.write('accession\tgenus\torganism\tassembly_level\tsource\n')
    for acc,vals in sorted(selected.items()):
        f.write('\t'.join((acc,*vals))+'\n')
(meta/'selected_accessions.txt').write_text(''.join(a+'\n' for a in sorted(selected)))
with (meta/'selection_summary.tsv').open('w') as f:
    f.write('genus\tmetadata_records\tselected\n')
    for row in stats: f.write('\t'.join(map(str,row))+'\n')
print('GENUS\tMETADATA\tSELECTED')
for row in stats: print(*row,sep='\t')
print(f'\nTOTAL ACCESSIONS (including existing controls): {len(selected)}')
print(f'ACCESSIONS: {meta}/selected_accessions.txt')
PY
  echo 'PLAN FINISHED. Inspect the accession count before downloading.'

elif [[ "$MODE" == download ]]; then
  require datasets; require python3; require timeout; require unzip
  [[ -s "$META/selected_accessions.txt" ]] || { echo 'Run plan first.' >&2; exit 1; }
  # Remove old batch *lists*; completed batch markers are validated against content.
  rm -f "$BATCH"/batch_*.txt
  python3 - "$META/selected_accessions.txt" "$BATCH" "$BATCH_SIZE" <<'PY'
import pathlib,sys
src=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2]); n=int(sys.argv[3])
acc=[x.strip() for x in src.read_text().splitlines() if x.strip()]
for i in range(0,len(acc),n):
    (out/f'batch_{i//n+1:03d}.txt').write_text('\n'.join(acc[i:i+n])+'\n')
print(f'Prepared {(len(acc)+n-1)//n} download batches for {len(acc)} accessions')
PY
  for LIST in "$BATCH"/batch_*.txt; do
    [[ -e "$LIST" ]] || continue
    STEM="${LIST%.txt}"; ZIP="${STEM}.zip"; DONE="${STEM}.done"
    TMPZIP="${STEM}.partial.zip"
    if [[ -e "$DONE" && -e "${STEM}.completed_accessions.txt" ]] \
        && cmp -s "$LIST" "${STEM}.completed_accessions.txt"; then
      echo "SKIP $(basename "$STEM") (complete)"; continue
    fi
    rm -f "$DONE"
    echo "DOWNLOAD $(basename "$STEM") ($(wc -l < "$LIST") genomes)"
    for TRY in 1 2; do
      rm -f "$TMPZIP"
      if timeout 1200s datasets download genome accession \
          --inputfile "$LIST" --include protein,gff3 \
          --filename "$TMPZIP" >"${STEM}.stdout" 2>"${STEM}.stderr" \
          && unzip -tq "$TMPZIP" >/dev/null; then
        mv "$TMPZIP" "$ZIP"
        unzip -oq "$ZIP" -d "$DATA"
        cp "$LIST" "${STEM}.completed_accessions.txt"
        touch "$DONE"
        break
      fi
      echo "  FAILED attempt $TRY; see ${STEM}.stderr" >&2
      sleep 5
    done
    [[ -e "$DONE" ]] || echo "WARNING: batch incomplete: $LIST" >&2
    sleep 1
  done
  echo "DOWNLOAD FINISHED; protein FASTAs: $(find "$DATA/ncbi_dataset/data" -name protein.faa 2>/dev/null | wc -l)"
  echo "Inspect incomplete batches with: ls $BATCH/batch_*.stderr"

elif [[ "$MODE" == analyze ]]; then
  require python3; require diamond
  [[ -s "$REF" ]] || { echo "Missing 13-gene reference: $REF" >&2; exit 1; }
  [[ -d "$DATA/ncbi_dataset/data" ]] || { echo 'No downloaded genomes; run download.' >&2; exit 1; }
  echo "== Building uniquely identified combined proteome =="
  python3 - "$DATA/ncbi_dataset/data" "$ANALYSIS" <<'PY'
import pathlib,sys
root=pathlib.Path(sys.argv[1]); out=pathlib.Path(sys.argv[2]); out.mkdir(parents=True,exist_ok=True)
fa=out/'combined_proteins.faa'; manifest=out/'downloaded_genomes.tsv'
proteins=0; assemblies=0
with fa.open('w') as dst,manifest.open('w') as mf:
    mf.write('accession\tprotein_fasta\tgff3\tproteins\n')
    for d in sorted(root.iterdir()):
        if not d.is_dir() or not d.name.startswith(('GCF_','GCA_')): continue
        src=d/'protein.faa'; gff=d/'genomic.gff'
        if not src.is_file() or src.stat().st_size==0: continue
        count=0
        with src.open(errors='replace') as inp:
            for line in inp:
                if line.startswith('>'):
                    pid=line[1:].split()[0]
                    dst.write(f'>{d.name}|{pid}\n')
                    count+=1
                else: dst.write(line)
        mf.write(f'{d.name}\t{src}\t{gff if gff.exists() else ""}\t{count}\n')
        proteins+=count; assemblies+=1
print(f'Combined {proteins:,} proteins from {assemblies} assemblies')
if assemblies==0: sys.exit('No protein FASTAs found')
PY
  echo '== Building DIAMOND database =='
  diamond makedb --in "$ANALYSIS/combined_proteins.faa" --db "$ANALYSIS/sol_discovery_db" \
    >"$ANALYSIS/makedb.log" 2>&1
  echo "== Searching all 13 sol proteins with $THREADS threads =="
  diamond blastp \
    --db "$ANALYSIS/sol_discovery_db" \
    --query "$REF" \
    --out "$ANALYSIS/diamond_hits.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \
    --sensitive --id 25 --query-cover 30 --evalue 1e-8 \
    --max-target-seqs 10000 --threads "$THREADS" \
    >"$ANALYSIS/diamond.log" 2>&1
  echo '== Mapping hits to GFF3 and identifying linked loci =='
  python3 - "$ANALYSIS" <<'PY'
import collections,csv,pathlib,re,sys,urllib.parse
out=pathlib.Path(sys.argv[1]); core={'solA','solB','solC','solD','solE','solF','solG','solH','solI','solK','solL'}
manifest={}
with (out/'downloaded_genomes.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'): manifest[r['accession']]=r
hits=collections.defaultdict(list); wanted=collections.defaultdict(set); raw_count=0
with (out/'diamond_hits.tsv').open() as f:
    for row in csv.reader(f,delimiter='\t'):
        if len(row)<12 or '|' not in row[1]: continue
        m=re.search(r'(?i)sol([A-M])',row[0])
        if not m: continue
        gene='sol'+m.group(1).upper(); acc,pid=row[1].split('|',1)
        if acc not in manifest: continue
        try: ident=float(row[2]); aln=int(row[3]); qlen=int(row[4]); bits=float(row[11])
        except ValueError: continue
        item=(gene,pid,ident,aln,qlen,bits)
        hits[acc].append(item); wanted[acc].add(pid); raw_count+=1

mapped=[]; unmapped=[]
for acc,rows in hits.items():
    path=pathlib.Path(manifest[acc]['gff3']) if manifest[acc]['gff3'] else None
    coords=collections.defaultdict(list)
    if path and path.is_file():
        with path.open(errors='replace') as f:
            for line in f:
                if line.startswith('#'): continue
                cols=line.rstrip('\n').split('\t')
                if len(cols)!=9 or cols[2]!='CDS': continue
                attrs={}
                for part in cols[8].split(';'):
                    if '=' in part:
                        k,v=part.split('=',1); attrs[k]=urllib.parse.unquote(v)
                ids=[attrs.get(k,'') for k in ('protein_id','Name','ID')]
                ids += [x.split(':',1)[1] for x in attrs.get('Dbxref','').split(',') if x.startswith(('Genbank:','RefSeq:'))]
                for pid in ids:
                    if pid in wanted[acc]:
                        try: coords[pid].append((cols[0],int(cols[3]),int(cols[4]),cols[6]))
                        except ValueError: pass
    for gene,pid,ident,aln,qlen,bits in rows:
        positions=coords.get(pid)
        if not positions:
            unmapped.append((acc,gene,pid,ident,aln,qlen,bits))
        else:
            for contig,start,end,strand in set(positions):
                mapped.append((acc,contig,start,end,strand,gene,pid,ident,aln,qlen,bits))

with (out/'mapped_hits.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['accession','contig','start','end','strand','gene','protein_id','identity_pct','aligned_aa','query_aa','bitscore']);w.writerows(mapped)
with (out/'unmapped_hits.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['accession','gene','protein_id','identity_pct','aligned_aa','query_aa','bitscore']);w.writerows(unmapped)

by_contig=collections.defaultdict(list)
for h in mapped: by_contig[(h[0],h[1])].append(h)
results=[]; max_span=100000
for (acc,contig),items in by_contig.items():
    items.sort(key=lambda x:(x[2],x[3]))
    best=None
    for i,first in enumerate(items):
        window=[]
        for h in items[i:]:
            if h[2]-first[2]>max_span: break
            if h[3]-first[2]>max_span: continue
            window.append(h)
        if not window: continue
        # Assign each protein to its highest-scoring sol query to avoid
        # double-counting multidomain proteins as multiple cluster genes.
        by_protein={}
        for h in window:
            if h[6] not in by_protein or h[10]>by_protein[h[6]][10]: by_protein[h[6]]=h
        by_gene={}
        for h in by_protein.values():
            if h[5] not in by_gene or h[10]>by_gene[h[5]][10]: by_gene[h[5]]=h
        genes=set(by_gene); ncore=len(genes&core); ntotal=len(genes)
        proteins=len({h[6] for h in by_gene.values()})
        start=min(h[2] for h in by_gene.values()); end=max(h[3] for h in by_gene.values())
        score=(ncore,ntotal,proteins,sum(h[10] for h in by_gene.values()),-(end-start))
        if best is None or score>best[0]: best=(score,start,end,genes,proteins,by_gene)
    if best:
        score,start,end,genes,proteins,by_gene=best
        ncore=len(genes&core)
        status=('HIGH_PRIORITY' if ncore>=9 and proteins>=8 else
                'FOLLOW_UP' if ncore>=7 and proteins>=6 else
                'POSSIBLE' if ncore>=5 and proteins>=5 else 'LOW')
        results.append((acc,contig,start,end,ncore,len(genes),proteins,','.join(sorted(genes&core)),','.join(sorted(genes)),status,round(score[3],1)))
results.sort(key=lambda x:(-x[4],-x[5],-x[6],-x[10],x[0]))
with (out/'linked_locus_candidates.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['accession','contig','window_start','window_end','shared_11_genes','all_13_genes','distinct_proteins','shared_genes','all_genes','priority','sum_best_bitscores']);w.writerows(results)
print(f'DIAMOND rows: {raw_count:,}; mapped rows: {len(mapped):,}; unmapped rows: {len(unmapped):,}')
print(f'Candidate contigs: {len(results):,}')
print('Priority counts:',dict(collections.Counter(x[9] for x in results)))
print('Top linked loci:')
for r in results[:15]: print('\t'.join(map(str,r[:10])))
print(f'RESULTS: {out}/linked_locus_candidates.tsv')
print('IMPORTANT: These are homology/synteny SCREENING calls, not validated sol cluster annotations.')
PY
  echo 'ANALYSIS FINISHED'
fi
