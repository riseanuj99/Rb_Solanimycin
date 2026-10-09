#!/bin/bash
#SBATCH --job-name=sol11_AU
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/sol11_AU_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/sol11_AU_%j.err

set -euo pipefail
export PATH="/scratch/al98750/Roux/envs/sol_hgt/bin:$PATH"
export ROOT=/scratch/al98750/Roux/07_sol_HGT
export OUT="$ROOT/07_phylogeny/sol11"
export THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "$OUT/alignments" "$OUT/targeted_AU" /scratch/al98750/Roux/logs
cd "$OUT"

# 1. Extract the 11 shared sol proteins, align, trim, concatenate.
python - <<'PY'
import csv, os, subprocess
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

root=Path(os.environ['ROOT']); out=Path(os.environ['OUT'])
old=root/'07_phylogeny'; aln=out/'alignments'
taxa=(old/'sol_tree_taxa.txt').read_text().splitlines()
genes=['solA','solB','solC','solD','solE','solF','solG','solH','solI','solK','solL']
assert len(taxa)==23 and len(set(taxa))==23, 'Expected 23 unique original taxa'

mapping=defaultdict(dict); low_coverage=[]
with (root/'05_results'/'sol_gene_hits.tsv').open() as f:
    for row in csv.DictReader(f,delimiter='\t'):
        a,g=row['accession'],row['gene']
        if a not in taxa or g not in genes: continue
        if g in mapping[a]: raise RuntimeError(f'Duplicate {a} {g}')
        mapping[a][g]=row['protein']
        if float(row['coverage'])<80:
            low_coverage.append((a,g,row['identity'],row['coverage'],row['protein']))
with (out/'qc_low_coverage.tsv').open('w') as f:
    f.write('accession\tgene\tidentity_pct\tcoverage_pct\tprotein\n')
    for row in low_coverage: f.write('\t'.join(row)+'\n')
print(f'QC warning: {len(low_coverage)} protein assignments have BLAST query coverage <80%; review before publication',flush=True)

seqs=defaultdict(dict)
for a in taxa:
    if set(mapping[a])!=set(genes):
        raise RuntimeError(f'Missing gene assignment in {a}: {set(genes)-set(mapping[a])}')
    wanted={pid:g for g,pid in mapping[a].items()}
    if len(wanted)!=len(genes): raise RuntimeError(f'Protein reused for different sol genes: {a}')
    faa=root/'03_proteomes'/f'{a}.faa'
    for record in SeqIO.parse(faa,'fasta'):
        if record.id in wanted:
            g=wanted[record.id]
            if g in seqs[a]: raise RuntimeError(f'Duplicate protein ID in {faa}: {record.id}')
            seqs[a][g]=str(record.seq)
    if set(seqs[a])!=set(genes):
        raise RuntimeError(f'Protein IDs not found in {faa}: {set(genes)-set(seqs[a])}')

concat={a:'' for a in taxa}; partitions=[]; pos=1
for g in genes:
    raw=aln/f'{g}.fa'; aligned=aln/f'{g}.aligned.fa'; trimmed=aln/f'{g}.trimmed.fa'
    with raw.open('w') as f:
        for a in taxa: f.write(f'>{a}\n{seqs[a][g]}\n')
    if not aligned.exists() or aligned.stat().st_size==0:
        temp=Path(str(aligned)+'.partial')
        with temp.open('w') as f:
            subprocess.run(['mafft','--auto',str(raw)],stdout=f,check=True)
        temp.replace(aligned)
    if not trimmed.exists() or trimmed.stat().st_size==0:
        temp=Path(str(trimmed)+'.partial')
        with temp.open('w') as f:
            subprocess.run(['trimal','-in',str(aligned),'-automated1'],stdout=f,check=True)
        temp.replace(trimmed)
    aligned_seqs=SeqIO.to_dict(SeqIO.parse(trimmed,'fasta'))
    if set(aligned_seqs)!=set(taxa): raise RuntimeError(f'Incorrect taxa in {g} alignment')
    lengths={len(r.seq) for r in aligned_seqs.values()}
    if len(lengths)!=1: raise RuntimeError(f'Unequal aligned sequence lengths for {g}')
    length=lengths.pop()
    if length<20: raise RuntimeError(f'{g} too short after trimming: {length}')
    for a in taxa: concat[a]+=str(aligned_seqs[a].seq)
    partitions.append((g,pos,pos+length-1)); pos+=length

with (out/'sol11_concat.fa').open('w') as f:
    for a in taxa: f.write(f'>{a}\n{concat[a]}\n')
with (out/'sol11_partitions.tsv').open('w') as f:
    f.write('gene\tstart\tend\n')
    for g,s,e in partitions: f.write(f'{g}\t{s}\t{e}\n')
with (out/'sol11_partitions.nex').open('w') as f:
    f.write('#nexus\nbegin sets;\n')
    for g,s,e in partitions: f.write(f'  charset {g} = {s}-{e};\n')
    f.write('end;\n')
(out/'sol11_taxa.txt').write_text('\n'.join(taxa)+'\n')
print(f'QC PASS: {len(taxa)} taxa x {len(genes)} genes = 253 proteins; {pos-1} aligned aa sites',flush=True)
print('Per-gene alignment lengths:',partitions,flush=True)
PY

IQ=iqtree3
ALN="$OUT/sol11_concat.fa"
MODEL=LG+F+G4

# 2. Unconstrained 11-gene ML tree with 1000 ultrafast bootstrap replicates.
if [[ ! -s "$OUT/sol11_ML.treefile" || ! -s "$OUT/sol11_ML.iqtree" ]]; then
    "$IQ" -s "$ALN" -m "$MODEL" -B 1000 -T "$THREADS" \
      --prefix "$OUT/sol11_ML"
fi

# 3. Build constraints based on organism names, not positional assumptions.
python - <<'PY'
import csv, os
from pathlib import Path
root=Path(os.environ['ROOT']); out=Path(os.environ['OUT']); d=out/'targeted_AU'
taxa=(out/'sol11_taxa.txt').read_text().splitlines()
with (root/'05_results'/'sol_loci.tsv').open() as f:
    org={r['accession']:r['organism'] for r in csv.DictReader(f,delimiter='\t')}
rb=[a for a in taxa if org[a]=='Rouxiella badensis']
rc=[a for a in taxa if org[a]=='Rouxiella chamberiensis']
ds=[a for a in taxa if org[a].startswith('Dickeya solani')]
pa=[a for a in taxa if org[a].startswith('Pantoea sp.')]
er=[a for a in taxa if org[a].startswith('Erwinia sp.')]
assert (len(rb),len(rc),len(ds),len(pa),len(er))==(18,2,1,1,1), 'Unexpected taxon groups'
def clade(x): return '('+','.join(x)+')'
# Unresolved within-group branches are intentionally free to optimize.
# The top-level split is the only alternative relationship of interest.
constraints={
 'RB_RC':f'(({clade(rb)},{clade(rc)}),({ds[0]},{pa[0]},{er[0]}));',
 'RB_PA':f'(({clade(rb)},{pa[0]}),({clade(rc)},{ds[0]},{er[0]}));'
}
for name,newick in constraints.items():
    (d/f'{name}.constraint.nwk').write_text(newick+'\n')
print('Constraints: RB+RC (species-compatible) versus RB+PA (alternative)',flush=True)
PY

# 4. Optimize each constrained topology on exactly the same alignment/model.
for H in RB_RC RB_PA; do
    PREFIX="$OUT/targeted_AU/$H"
    if [[ ! -s "$PREFIX.treefile" || ! -s "$PREFIX.iqtree" ]]; then
        "$IQ" -s "$ALN" -m "$MODEL" \
          -g "$OUT/targeted_AU/$H.constraint.nwk" \
          -T "$THREADS" --prefix "$PREFIX"
    fi
done

# 5. Validate the focal unrooted splits and remove any duplicate topologies.
python - <<'PY'
import csv, os
from pathlib import Path
from Bio import Phylo
root=Path(os.environ['ROOT']); out=Path(os.environ['OUT']); d=out/'targeted_AU'
taxa=(out/'sol11_taxa.txt').read_text().splitlines(); alltaxa=set(taxa)
with (root/'05_results'/'sol_loci.tsv').open() as f:
    org={r['accession']:r['organism'] for r in csv.DictReader(f,delimiter='\t')}
rb={a for a in taxa if org[a]=='Rouxiella badensis'}
rc={a for a in taxa if org[a]=='Rouxiella chamberiensis'}
pa={a for a in taxa if org[a].startswith('Pantoea sp.')}
def bipartitions(tree):
    splits=set()
    for clade in tree.find_clades():
        side=frozenset(t.name for t in clade.get_terminals())
        other=frozenset(alltaxa-side)
        if len(side)>1 and len(other)>1:
            splits.add(frozenset((side,other)))
    return frozenset(splits)
def has_split(tree, group):
    return frozenset((frozenset(group),frozenset(alltaxa-group))) in bipartitions(tree)
files=[('unconstrained',out/'sol11_ML.treefile'),
       ('RB_RC',d/'RB_RC.treefile'),('RB_PA',d/'RB_PA.treefile')]
seen={}; unique=[]; table=[]
for label,path in files:
    tree=Phylo.read(path,'newick')
    found={t.name for t in tree.get_terminals()}
    if found!=alltaxa: raise RuntimeError(f'{label}: incorrect taxa: {found^alltaxa}')
    if label=='RB_RC' and not has_split(tree,rb|rc):
        raise RuntimeError('RB_RC constraint NOT recovered; stop rather than report invalid AU test')
    if label=='RB_PA' and not has_split(tree,rb|pa):
        raise RuntimeError('RB_PA constraint NOT recovered; stop rather than report invalid AU test')
    sig=bipartitions(tree)
    if sig in seen:
        treeid=seen[sig]; status='duplicate_topology'
    else:
        treeid=len(unique)+1; seen[sig]=treeid; unique.append(path.read_text().strip()); status='unique'
    table.append((label,treeid,status,str(path)))
if len(unique)<2: raise RuntimeError('Fewer than two distinct candidate topologies')
(d/'candidate_trees.nwk').write_text('\n'.join(unique)+'\n')
with (d/'candidate_order.tsv').open('w') as f:
    f.write('hypothesis\tAU_TreeID\tstatus\toptimized_tree\n')
    for row in table: f.write('\t'.join(map(str,row))+'\n')
print('Targeted AU tree order:',table,flush=True)
PY

# 6. Targeted AU test with 10,000 RELL bootstrap replicates.
if [[ ! -s "$OUT/targeted_AU/focal_AU_10k.iqtree" ]]; then
    "$IQ" -s "$ALN" -m "$MODEL" \
      -z "$OUT/targeted_AU/candidate_trees.nwk" \
      -n 0 -zb 10000 -au -T "$THREADS" \
      --prefix "$OUT/targeted_AU/focal_AU_10k"
fi

# 7. Individual-gene trees; skip those already completed on reruns.
for G in solA solB solC solD solE solF solG solH solI solK solL; do
    PREFIX="$OUT/alignments/${G}_ML"
    if [[ ! -s "$PREFIX.treefile" || ! -s "$PREFIX.iqtree" ]]; then
        "$IQ" -s "$OUT/alignments/${G}.trimmed.fa" \
          -m "$MODEL" -B 1000 -T "$THREADS" --prefix "$PREFIX"
    fi
done

printf '\n========== COMPLETE: 11-GENE PHYLOGENY + TARGETED AU =========='\nprintf '\nAlignment: %s\nTree: %s\n' "$ALN" "$OUT/sol11_ML.treefile"
printf 'Candidate order: %s\nAU report: %s\n' \
  "$OUT/targeted_AU/candidate_order.tsv" "$OUT/targeted_AU/focal_AU_10k.iqtree"
echo 'NOTE: AU p-values test candidate topologies, not the probability of HGT.'
