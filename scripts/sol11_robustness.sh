#!/bin/bash
#SBATCH --job-name=sol11_robust
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/al98750/Roux/logs/sol11_robust_%j.out
#SBATCH --error=/scratch/al98750/Roux/logs/sol11_robust_%j.err

# Robustness analyses for the 23-taxon sol HGT dataset.
# Reads previous sol11 results; never overwrites the original trees or AU test.
# 1. Summarize 11 existing individual-gene trees with root-independent splits.
# 2. Reconstruct 11-gene phylogeny with 11 proportional partitions.
# 3. Reconstruct 10-gene phylogeny excluding solA, with 10 partitions.
# 4. Run focal RB+RC vs RB+Pantoea AU tests (10,000 RELL replicates)
#    under each partitioned model, after constrained ML optimization.

set -euo pipefail
export PATH="/scratch/al98750/Roux/envs/sol_hgt/bin:$PATH"
export ROOT="/scratch/al98750/Roux/07_sol_HGT"
export BASE="$ROOT/07_phylogeny/sol11"
export OUT="$BASE/robustness"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
IQ="/scratch/al98750/Roux/envs/sol_hgt/bin/iqtree3"

mkdir -p "$OUT" "$OUT/partitioned11" "$OUT/noSolA_partitioned10" /scratch/al98750/Roux/logs

command -v python >/dev/null || { echo 'Python not found'; exit 1; }
[[ -x "$IQ" ]] || { echo "IQ-TREE executable not found: $IQ"; exit 1; }

echo "========== 1. INPUT QC, GENE-TREE CONCORDANCE, NO-solA ALIGNMENT =========="
python - <<'PY'
import csv
import os
from collections import Counter
from pathlib import Path
from Bio import SeqIO, Phylo

root = Path(os.environ['ROOT'])
base = Path(os.environ['BASE'])
out = Path(os.environ['OUT'])
allgenes = ['solA','solB','solC','solD','solE','solF','solG','solH','solI','solK','solL']
keptgenes = [g for g in allgenes if g != 'solA']

taxa = (base/'sol11_taxa.txt').read_text().splitlines()
if len(taxa) != 23 or len(set(taxa)) != 23:
    raise RuntimeError('Expected 23 distinct taxa from completed sol11 analysis')
alltaxa = frozenset(taxa)

with (root/'05_results/sol_loci.tsv').open() as fh:
    org = {r['accession']:r['organism'] for r in csv.DictReader(fh, delimiter='\t')}
rb = frozenset(a for a in taxa if org[a] == 'Rouxiella badensis')
rc = frozenset(a for a in taxa if org[a] == 'Rouxiella chamberiensis')
ds = frozenset(a for a in taxa if org[a].startswith('Dickeya solani'))
pa = frozenset(a for a in taxa if org[a].startswith('Pantoea sp.'))
er = frozenset(a for a in taxa if org[a].startswith('Erwinia sp.'))
if tuple(map(len, (rb,rc,ds,pa,er))) != (18,2,1,1,1):
    raise RuntimeError('Unexpected taxon composition; verify original phylogeny')

def splits_with_support(tree):
    """Root-independent bipartitions, with edge supports when available."""
    tips = frozenset(t.name for t in tree.get_terminals())
    if tips != alltaxa:
        raise RuntimeError(f'Tree taxa mismatch: {sorted(tips ^ alltaxa)}')
    found = {}
    for node in tree.find_clades():
        side = frozenset(t.name for t in node.get_terminals())
        complement = alltaxa - side
        if len(side) <= 1 or len(complement) <= 1:
            continue
        split = frozenset((side, complement))
        support = node.confidence
        if split not in found or (found[split] is None and support is not None):
            found[split] = support
    return found

def support_for(splits, group):
    target = frozenset((frozenset(group), alltaxa - frozenset(group)))
    if target not in splits:
        return 'absent'
    val = splits[target]
    return 'present_no_support' if val is None else f'{float(val):g}'

concat = {a: '' for a in taxa}
partitions = []
pos = 1
rows = []
for gene in allgenes:
    fasta = base/'alignments'/f'{gene}.trimmed.fa'
    treefile = base/'alignments'/f'{gene}_ML.treefile'
    if not fasta.is_file() or not treefile.is_file():
        raise FileNotFoundError(f'Missing completed input for {gene}: {fasta} or {treefile}')
    records = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    if set(records) != set(taxa):
        raise RuntimeError(f'{gene}: expected 23 taxa in alignment')
    lengths = {len(rec.seq) for rec in records.values()}
    if len(lengths) != 1:
        raise RuntimeError(f'{gene}: unequal alignment lengths')
    length = lengths.pop()
    if length < 20:
        raise RuntimeError(f'{gene}: alignment unexpectedly short')
    tree = Phylo.read(treefile, 'newick')
    splits = splits_with_support(tree)
    rbpa = support_for(splits, rb|pa)
    rbrc = support_for(splits, rb|rc)
    if rbpa != 'absent' and rbrc != 'absent':
        raise RuntimeError(f'{gene}: incompatible RB+PA and RB+RC splits both detected')
    rows.append({
        'gene': gene, 'alignment_aa': length,
        'unique_sequences': len(set(str(records[a].seq) for a in taxa)),
        'RB_monophyly_support': support_for(splits, rb),
        'RC_monophyly_support': support_for(splits, rc),
        'RB_PA_split_support': rbpa,
        'RB_RC_split_support': rbrc,
        'DS_ER_split_support': support_for(splits, ds|er),
        'focal_grouping': ('RB_PA' if rbpa != 'absent' else
                           'RB_RC' if rbrc != 'absent' else 'neither')
    })
    if gene in keptgenes:
        for a in taxa:
            concat[a] += str(records[a].seq)
        partitions.append((gene, pos, pos+length-1))
        pos += length

with (out/'gene_tree_focal_splits.tsv').open('w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)

# 10-gene concatenation uses the EXACT same trimmed gene alignments as the 11-gene run.
aln10 = out/'noSolA_partitioned10'/'sol10_concat.fa'
with aln10.open('w') as fh:
    for a in taxa:
        fh.write(f'>{a}\n{concat[a]}\n')

part10 = out/'noSolA_partitioned10'/'sol10_partitions.nex'
with part10.open('w') as fh:
    fh.write('#nexus\nbegin sets;\n')
    for gene, start, end in partitions:
        fh.write(f'  charset {gene} = {start}-{end};\n')
    fh.write('end;\n')
with (out/'noSolA_partitioned10'/'sol10_partitions.tsv').open('w') as fh:
    fh.write('gene\tstart\tend\n')
    for gene, start, end in partitions:
        fh.write(f'{gene}\t{start}\t{end}\n')

full_records = SeqIO.to_dict(SeqIO.parse(base/'sol11_concat.fa', 'fasta'))
if set(full_records) != set(taxa):
    raise RuntimeError('Original concatenated 11-gene alignment has unexpected taxa')
if any(len(rec.seq) != 12987 for rec in full_records.values()):
    raise RuntimeError('Original 11-gene alignment length changed; inspect before running')
if len(next(iter(concat.values()))) != 7840:
    raise RuntimeError('Expected 7,840 aa after excluding solA')

print('Original: 23 taxa, 11 genes, 12,987 aa', flush=True)
print('Without solA: 23 taxa, 10 genes, 7,840 aa', flush=True)
print('Individual-gene focal grouping counts:',
      dict(Counter(r['focal_grouping'] for r in rows)), flush=True)
print('Gene-tree table:', out/'gene_tree_focal_splits.tsv', flush=True)
PY

# Edge-linked proportional partition model: -p allows each gene its own
# relative evolutionary rate and substitution-model parameters.
# LG+F+G4 is fixed across partitions for a controlled sensitivity analysis.
# ModelFinder-per-gene can be added later as a separate model-selection check.
MODEL="LG+F+G4"
SEED=20261012

run_partitioned() {
    local dir="$1" aln="$2" part="$3" prefix="$4"
    if [[ -s "$prefix.treefile" && -s "$prefix.iqtree" ]]; then
        echo "SKIP completed tree: $prefix"
    else
        echo "INFER partitioned tree: $prefix"
        "$IQ" -s "$aln" -p "$part" -m "$MODEL" \
            -B 1000 -T "$THREADS" -seed "$SEED" --prefix "$prefix"
    fi
}

run_focal_au() {
    local dir="$1" aln="$2" part="$3" best="$4"
    local h prefix
    for h in RB_RC RB_PA; do
        prefix="$dir/targeted_AU/$h"
        if [[ -s "$prefix.treefile" && -s "$prefix.iqtree" ]]; then
            echo "SKIP completed constrained tree: $prefix"
        else
            echo "INFER constrained tree: $h ($dir)"
            "$IQ" -s "$aln" -p "$part" -m "$MODEL" \
                -g "$BASE/targeted_AU/$h.constraint.nwk" \
                -T "$THREADS" -seed "$SEED" --prefix "$prefix"
        fi
    done

    # Confirm BOTH focal constraints were honored, regardless of tree rooting;
    # deduplicate tree topologies before computing the AU confidence set.
    RUN_DIR="$dir" BEST_TREE="$best" python - <<'PY'
import csv
import os
from pathlib import Path
from Bio import Phylo

root = Path(os.environ['ROOT'])
base = Path(os.environ['BASE'])
dir = Path(os.environ['RUN_DIR'])
d = dir/'targeted_AU'
taxa = (base/'sol11_taxa.txt').read_text().splitlines()
alltaxa = frozenset(taxa)
with (root/'05_results/sol_loci.tsv').open() as fh:
    org = {r['accession']:r['organism'] for r in csv.DictReader(fh,delimiter='\t')}
rb = frozenset(a for a in taxa if org[a]=='Rouxiella badensis')
rc = frozenset(a for a in taxa if org[a]=='Rouxiella chamberiensis')
pa = frozenset(a for a in taxa if org[a].startswith('Pantoea sp.'))

def splits(tree):
    if frozenset(t.name for t in tree.get_terminals()) != alltaxa:
        raise RuntimeError('Incorrect taxon set in ' + str(dir))
    result=set()
    for node in tree.find_clades():
        side=frozenset(t.name for t in node.get_terminals())
        other=alltaxa-side
        if len(side)>1 and len(other)>1:
            result.add(frozenset((side,other)))
    return frozenset(result)

def has_group(s,group):
    return frozenset((frozenset(group), alltaxa-frozenset(group))) in s

candidates=[('unconstrained',Path(os.environ['BEST_TREE'])),
            ('RB_RC',d/'RB_RC.treefile'),
            ('RB_PA',d/'RB_PA.treefile')]
seen={}; unique=[]; table=[]
for name,path in candidates:
    tree=Phylo.read(path,'newick')
    sig=splits(tree)
    if name=='RB_RC' and not has_group(sig,rb|rc):
        raise RuntimeError('RB_RC constrained tree lacks required unrooted split')
    if name=='RB_PA' and not has_group(sig,rb|pa):
        raise RuntimeError('RB_PA constrained tree lacks required unrooted split')
    if sig in seen:
        idx=seen[sig]; status='duplicate_topology'
    else:
        idx=len(unique)+1; seen[sig]=idx
        unique.append(path.read_text().strip()); status='unique'
    table.append((name,idx,status,str(path)))
if len(unique)<2:
    raise RuntimeError('Fewer than two distinct topologies; AU test cannot proceed')
(d/'candidate_trees.nwk').write_text('\n'.join(unique)+'\n')
with (d/'candidate_order.tsv').open('w') as fh:
    fh.write('hypothesis\tAU_TreeID\tstatus\toptimized_tree\n')
    for r in table: fh.write('\t'.join(map(str,r))+'\n')
print('Validated focal hypotheses:',table,flush=True)
PY

    if [[ -s "$dir/targeted_AU/focal_AU_10k.iqtree" ]]; then
        echo "SKIP completed AU: $dir"
    else
        echo "RUN 10,000-replicate partitioned AU: $dir"
        "$IQ" -s "$aln" -p "$part" -m "$MODEL" \
            -z "$dir/targeted_AU/candidate_trees.nwk" \
            -n 0 -zb 10000 -au -T "$THREADS" -seed "$SEED" \
            --prefix "$dir/targeted_AU/focal_AU_10k"
    fi
}

DIR11="$OUT/partitioned11"
DIR10="$OUT/noSolA_partitioned10"
mkdir -p "$DIR11/targeted_AU" "$DIR10/targeted_AU"

ALN11="$BASE/sol11_concat.fa"
PART11="$BASE/sol11_partitions.nex"
ALN10="$DIR10/sol10_concat.fa"
PART10="$DIR10/sol10_partitions.nex"

for f in "$ALN11" "$PART11" "$ALN10" "$PART10" \
    "$BASE/targeted_AU/RB_RC.constraint.nwk" \
    "$BASE/targeted_AU/RB_PA.constraint.nwk"; do
    [[ -s "$f" ]] || { echo "Required input missing: $f"; exit 1; }
done

echo "========== 2. PARTITIONED 11-GENE PHYLOGENY =========="
run_partitioned "$DIR11" "$ALN11" "$PART11" "$DIR11/sol11_partitioned_ML"
echo "========== 3. TARGETED AU: PARTITIONED 11-GENE =========="
run_focal_au "$DIR11" "$ALN11" "$PART11" "$DIR11/sol11_partitioned_ML.treefile"

echo "========== 4. PARTITIONED 10-GENE PHYLOGENY (EXCLUDE solA) =========="
run_partitioned "$DIR10" "$ALN10" "$PART10" "$DIR10/sol10_partitioned_ML"
echo "========== 5. TARGETED AU: 10-GENE NO-solA =========="
run_focal_au "$DIR10" "$ALN10" "$PART10" "$DIR10/sol10_partitioned_ML.treefile"

echo "========== 6. SUMMARIZE TREE SPLITS AND AU TESTS =========="
python - <<'PY'
import csv
import os
import re
from pathlib import Path
from Bio import Phylo

root=Path(os.environ['ROOT']); base=Path(os.environ['BASE']); out=Path(os.environ['OUT'])
taxa=(base/'sol11_taxa.txt').read_text().splitlines()
alltaxa=frozenset(taxa)
with (root/'05_results/sol_loci.tsv').open() as fh:
    org={r['accession']:r['organism'] for r in csv.DictReader(fh,delimiter='\t')}
rb=frozenset(a for a in taxa if org[a]=='Rouxiella badensis')
rc=frozenset(a for a in taxa if org[a]=='Rouxiella chamberiensis')
pa=frozenset(a for a in taxa if org[a].startswith('Pantoea sp.'))

def support(tree,group):
    target=frozenset((frozenset(group),alltaxa-frozenset(group)))
    values=[]
    for node in tree.find_clades():
        side=frozenset(t.name for t in node.get_terminals())
        other=alltaxa-side
        if len(side)>1 and len(other)>1 and frozenset((side,other))==target:
            values.append(node.confidence)
    if not values: return 'absent'
    return str(next((float(v) for v in values if v is not None),'present_no_support'))

def au_pvalues(report):
    lines=report.read_text().splitlines()
    try: start=next(i for i,l in enumerate(lines) if l.strip()=='USER TREES')
    except StopIteration: raise RuntimeError('USER TREES section absent in '+str(report))
    pvals={}
    for line in lines[start:]:
        if line.strip().startswith('deltaL') and pvals: break
        m=re.match(r'^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+',line)
        if not m: continue
        tokens=line.split()
        if tokens[-1] not in ('+','-'):
            raise RuntimeError('Unexpected AU table format: '+line)
        pvals[int(tokens[0])]=float(tokens[-2])
    if not pvals: raise RuntimeError('No AU p-values found in '+str(report))
    return pvals

rows=[]
for name,subdir,treefile,sites in [
    ('partitioned_11','partitioned11','sol11_partitioned_ML.treefile',12987),
    ('partitioned_no_solA_10','noSolA_partitioned10','sol10_partitioned_ML.treefile',7840),
]:
    d=out/subdir
    tree=Phylo.read(d/treefile,'newick')
    if frozenset(t.name for t in tree.get_terminals())!=alltaxa:
        raise RuntimeError('Final tree has wrong taxon set: '+name)
    pvals=au_pvalues(d/'targeted_AU/focal_AU_10k.iqtree')
    with (d/'targeted_AU/candidate_order.tsv').open() as fh:
        ids={r['hypothesis']:int(r['AU_TreeID']) for r in csv.DictReader(fh,delimiter='\t')}
    rows.append({
        'analysis':name,'taxa':23,'aligned_aa':sites,
        'RB_PA_bootstrap':support(tree,rb|pa),
        'RB_RC_bootstrap':support(tree,rb|rc),
        'RB_bootstrap':support(tree,rb),
        'RC_bootstrap':support(tree,rc),
        'AU_unconstrained':f"{pvals[ids['unconstrained']]:.5g}",
        'AU_RB_RC':f"{pvals[ids['RB_RC']]:.5g}",
        'AU_RB_PA':f"{pvals[ids['RB_PA']]:.5g}",
    })
with (out/'robustness_summary.tsv').open('w',newline='') as fh:
    writer=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t')
    writer.writeheader(); writer.writerows(rows)
print('ROBUSTNESS SUMMARY:')
with (out/'robustness_summary.tsv').open() as fh:
    print(fh.read(),flush=True)
print('INDIVIDUAL GENE TREE SUMMARY:')
with (out/'gene_tree_focal_splits.tsv').open() as fh:
    print(fh.read(),flush=True)
PY

echo '========== COMPLETE: GENE-TREE CONCORDANCE + PARTITIONED AND NO-solA AU =========='
echo "Results: $OUT"
