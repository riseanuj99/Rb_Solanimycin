#!/usr/bin/env python3
"""Create publication-oriented figures and tables from existing Sapelo2 sol HGT outputs.

Dependencies: biopython, matplotlib (standard library otherwise).
Usage:
 python make_sol_hgt_figures.py \
   --root /scratch/al98750/Roux/07_sol_HGT/07_phylogeny \
   --out /scratch/al98750/Roux/07_sol_HGT/10_manuscript_figures

No alignments, trees, or phylogenetic analyses are modified or rerun.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

RB_ACCESSIONS = {
    'GCA_026967495', 'GCF_002093665', 'GCF_011067265',
    'GCF_014859005', 'GCF_020739645', 'GCF_020739655',
    'GCF_020740165', 'GCF_020740185', 'GCF_020740205',
    'GCF_020740225', 'GCF_020740305', 'GCF_026967515',
    'GCF_032083745', 'GCF_037153295', 'GCF_049580565',
    'GCF_049582605', 'GCF_049590195', 'GCF_963516055',
}
RC_ACCESSIONS = {'GCF_000951135', 'GCF_026967475'}
COLORS = {
    'RB': '#167b78',
    'RC': '#4875a9',
    'PA': '#b86d1d',
    'ER': '#7d61a8',
    'DS': '#555b64',
    'Unknown': '#9a9a9a',
}
GROUP_LABELS = {
    'RB': 'R. badensis', 'RC': 'R. chamberiensis',
    'PA': 'Pantoea sp. A4', 'ER': 'Erwinia sp. ErVv1',
    'DS': 'D. solani MK10', 'Unknown': 'Unclassified'
}
GROUP_ORDER = {'RB': 0, 'RC': 1, 'PA': 2, 'ER': 3, 'DS': 4, 'Unknown': 5}


def classify(name: str) -> str:
    match = re.search(r'(?:GCA|GCF)_\d+', name)
    acc = match.group(0) if match else ''
    if acc in RB_ACCESSIONS or 'badensis' in name.lower():
        return 'RB'
    if acc in RC_ACCESSIONS or 'chamberiensis' in name.lower():
        return 'RC'
    if acc == 'GCA_000365285' or 'dsolani' in name.lower() or 'dickeya' in name.lower():
        return 'DS'
    if acc == 'GCA_000295955' or 'pantoea' in name.lower():
        return 'PA'
    if acc == 'GCA_900068895' or 'erwinia' in name.lower():
        return 'ER'
    return 'Unknown'


def taxon_display(name: str) -> str:
    group = classify(name)
    acc = re.search(r'(?:GCA|GCF)_\d+(?:\.\d+)?', name)
    identifier = acc.group(0) if acc else name
    return f'{group}  {identifier}'


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(f'Required result not found: {path}')
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_table(path: Path, rows: list[dict], columns: list[str], sep: str):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter=sep, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def save_figure(fig, out: Path, stem: str):
    for ext in ['pdf', 'svg', 'png']:
        kwargs = {'dpi': 300} if ext == 'png' else {}
        fig.savefig(out / f'{stem}.{ext}', bbox_inches='tight', facecolor='white', **kwargs)
    plt.close(fig)


def au_results(path: Path) -> dict[str, float]:
    """Parse p-AU column using the IQ-TREE USER TREES table, not guessed positions."""
    text = path.read_text()
    section = text.split('USER TREES', 1)
    if len(section) != 2:
        raise ValueError(f'USER TREES section not found: {path}')
    table = section[1].split('ALISIM COMMAND', 1)[0]
    values = {}
    for line in table.splitlines():
        m = re.match(r'^\s*([123])\s+(-?\d+(?:\.\d+)?)\s+', line)
        if not m:
            continue
        tokens = line.split()
        # IQ-TREE appends a + or - confidence-set indicator after p-AU.
        p_token = tokens[-2] if tokens[-1] in {'+', '-'} else tokens[-1]
        try:
            values[int(tokens[0])] = float(p_token)
        except ValueError as exc:
            raise ValueError(f'Could not parse AU p-value from: {line}') from exc
    if set(values) != {1, 2, 3}:
        raise ValueError(f'Expected three AU candidates in {path}, found {values}')
    # Do not assume candidate order; use candidate_order.tsv if it exists.
    order_file = path.parent / 'candidate_order.tsv'
    if order_file.exists():
        mapping = {}
        for row in read_tsv(order_file):
            key = row.get('hypothesis', '').strip()
            idx = int(row.get('AU_TreeID', 0))
            mapping[key] = values[idx]
        if {'unconstrained', 'RB_RC', 'RB_PA'} <= set(mapping):
            return mapping
    raise ValueError(f'Could not verify AU tree-to-hypothesis mapping in {order_file}')


def leaves(tree) -> list[str]:
    return [clade.name for clade in tree.get_terminals()]


def validate_taxa(trees: dict[str, object]):
    baseline = None
    for label, tree in trees.items():
        names = leaves(tree)
        groups = Counter(classify(x) for x in names)
        if len(names) != 23 or len(set(names)) != 23:
            raise ValueError(f'{label}: expected 23 distinct terminal taxa; found {len(names)} tips / {len(set(names))} distinct')
        expected = {'RB': 18, 'RC': 2, 'PA': 1, 'ER': 1, 'DS': 1}
        if dict(groups) != expected:
            raise ValueError(f'{label}: unexpected taxon classification {dict(groups)}; expected {expected}.\n'
                             f'Please check taxon IDs:\n' + '\n'.join(names))
        if baseline is None:
            baseline = set(names)
        elif set(names) != baseline:
            raise ValueError(f'{label}: tip labels differ from species tree.\nMissing: {sorted(baseline-set(names))}\nExtra: {sorted(set(names)-baseline)}')
        print(f'  {label}: 23 taxa, groups {dict(groups)}')


def split_support(tree, focal: set[str]) -> float | None:
    """Find branch support for an UNROOTED bipartition; allow either split side."""
    all_names = set(leaves(tree))
    for clade in tree.find_clades(order='preorder'):
        if clade is tree.root:
            continue
        descendants = {x.name for x in clade.get_terminals()}
        if descendants == focal or all_names - descendants == focal:
            if clade.confidence is not None:
                return float(clade.confidence)
    return None


def focal_set(tree, *groups):
    return {name for name in leaves(tree) if classify(name) in groups}


def sort_clades(clade):
    for sub in clade.clades:
        sort_clades(sub)
    def key(child):
        terms = child.get_terminals()
        return (min(GROUP_ORDER[classify(x.name)] for x in terms),
                min(x.name for x in terms))
    clade.clades.sort(key=key)


def draw_cladogram(ax, source_tree, title, highlight_groups, support=None, font_size=7.1):
    tree = deepcopy(source_tree)
    ds = [x for x in tree.get_terminals() if classify(x.name) == 'DS']
    if len(ds) != 1:
        raise ValueError('Expected exactly one D. solani tip for display rooting')
    # This is an arbitrary VISUAL root. Do not interpret as an inferred outgroup.
    tree.root_with_outgroup(ds[0])
    sort_clades(tree.root)
    x = tree.depths(unit_branch_lengths=True)
    terminals = tree.get_terminals()
    y = {leaf: i for i, leaf in enumerate(terminals)}
    def fill_y(clade):
        if clade not in y:
            ys = [fill_y(c) for c in clade.clades]
            y[clade] = (min(ys) + max(ys)) / 2
        return y[clade]
    fill_y(tree.root)
    max_depth = max(x.values())
    target = focal_set(tree, *highlight_groups)
    highlight = None
    for clade in tree.find_clades(order='preorder'):
        if {t.name for t in clade.get_terminals()} == target:
            highlight = clade
            break
    for clade in tree.find_clades(order='preorder'):
        if not clade.clades:
            continue
        ys = [y[c] for c in clade.clades]
        ax.plot([x[clade], x[clade]], [min(ys), max(ys)], lw=0.9, color='#374151', zorder=2)
        for child in clade.clades:
            color = '#c25136' if child is highlight else '#374151'
            width = 2.3 if child is highlight else 0.9
            ax.plot([x[clade], x[child]], [y[child], y[child]], lw=width, color=color, zorder=3)
    for leaf in terminals:
        ax.text(x[leaf] + 0.12, y[leaf], taxon_display(leaf.name),
                fontsize=font_size, color=COLORS[classify(leaf.name)], va='center', clip_on=False)
    if highlight is not None and support is not None:
        ax.scatter([x[highlight]], [y[highlight]], s=30, color='#c25136', zorder=4)
    ax.set_xlim(-0.2, max_depth + 5.5)
    ax.set_ylim(len(terminals) - 0.4, -1.1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    suffix = f'Focal split UFBoot: {support:g}%' if support is not None else 'Focal split not recovered'
    ax.set_title(f'{title}\n{suffix}', loc='left', fontsize=12, fontweight='bold', pad=13)


def figure1(species_tree, sol_tree, out):
    species_support = split_support(species_tree, focal_set(species_tree, 'RB', 'RC'))
    sol_support = split_support(sol_tree, focal_set(sol_tree, 'RB', 'PA'))
    fig, axes = plt.subplots(1, 2, figsize=(21, 10.5), gridspec_kw={'wspace': 0.22})
    draw_cladogram(axes[0], species_tree, 'A   Core-genome species phylogeny (Rb + Rc)', ('RB', 'RC'), species_support)
    draw_cladogram(axes[1], sol_tree, 'B   Partitioned 11-gene sol phylogeny (Rb + Pa)', ('RB', 'PA'), sol_support)
    handles = [Line2D([0], [0], marker='o', markersize=8, linestyle='', color=COLORS[g], label=GROUP_LABELS[g])
               for g in ['RB', 'RC', 'PA', 'ER', 'DS']]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.045), ncol=5, frameon=False, fontsize=9)
    fig.text(0.5, 0.015,
             '23 genomes | Cladograms: branch lengths not to scale | Rooted on D. solani for visualization only (not an inferred outgroup).',
             ha='center', fontsize=9, color='#525b65')
    fig.subplots_adjust(bottom=0.12, top=0.89, left=0.03, right=0.97)
    save_figure(fig, out, 'Figure1_species_vs_partitioned_sol_tree')
    return species_support, sol_support


def figure2(gene_rows, out):
    labels = [r['gene'] for r in gene_rows]
    lengths = [int(r['alignment_aa']) for r in gene_rows]
    values = [float(r['RB_PA_split_support']) if r['RB_PA_split_support'] not in {'absent', 'NA', ''} else None for r in gene_rows]
    y = list(range(len(labels)))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 6.8), sharey=True,
                                   gridspec_kw={'width_ratios': [1.25, 0.9], 'wspace': 0.11})
    for i, value in enumerate(values):
        if value is None:
            ax1.barh(i, 0.6, color='#b5b9bf', height=0.67)
            ax1.text(2.5, i, 'Not recovered', fontsize=9, va='center', color='#666666')
        else:
            col = '#167b78' if value >= 80 else '#c28b37'
            ax1.barh(i, value, color=col, height=0.67)
            ax1.text(min(value + 2, 106), i, f'{value:g}%', fontsize=9, va='center')
    ax1.set_xlim(0, 119)
    ax1.set_yticks(y, labels=[f'$\\mathit{{{s}}}$' for s in labels])
    ax1.invert_yaxis()
    ax1.set_xlabel('UFBoot support for Rb + Pantoea (%)')
    ax1.set_title('A   Individual-gene topology support', loc='left', fontsize=12, fontweight='bold')
    ax1.axvline(80, color='#6b7280', lw=1, ls='--', alpha=0.55)
    ax1.set_xticks([0, 20, 40, 60, 80, 100])
    ax1.grid(axis='x', color='#e7eaee', lw=0.6)
    ax1.set_axisbelow(True)
    ax2.barh(y, lengths, height=0.67, color=['#b86d1d' if n == 'solA' else '#94a8b9' for n in labels])
    for i, length in enumerate(lengths):
        ax2.text(length + 55, i, f'{length:,}', va='center', fontsize=9)
    ax2.set_xlim(0, max(lengths) * 1.25)
    ax2.set_xlabel('Trimmed alignment length (aa)')
    ax2.set_title('B   Gene alignment contribution', loc='left', fontsize=12, fontweight='bold')
    ax2.grid(axis='x', color='#e7eaee', lw=0.6)
    ax2.set_axisbelow(True)
    for ax in (ax1, ax2):
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.tick_params(axis='y', length=0)
    fig.text(0.03, 0.025,
             '“Not recovered” is not 0% support. All genes have 23 taxa; linked genes are not independent HGT events.',
             fontsize=9, color='#4b5563')
    fig.subplots_adjust(left=0.115, right=0.975, bottom=0.14, top=0.87)
    save_figure(fig, out, 'Figure2_individual_gene_support')


def figure3(rows, out):
    labels = ['Unpartitioned\n11 genes', 'Partitioned\n11 genes', 'Partitioned\n10 genes (no solA)']
    boots = [float(r['RB_PA_bootstrap']) if r['RB_PA_bootstrap'] not in {'absent', 'NA', ''} else None for r in rows]
    pvals = [float(r['AU_RB_RC']) for r in rows]
    logp = [-math.log10(max(p, 1e-300)) for p in pvals]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.7, 5.5), gridspec_kw={'wspace': 0.32})
    xs = list(range(3))
    ax1.bar(xs, [v if v is not None else 0 for v in boots], color=['#8094aa', '#167b78', '#167b78'], width=0.58)
    if boots[0] is None:
        ax1.text(0, 12, 'Not yet verified', ha='center', fontsize=9, color='#555b64', rotation=90)
    ax1.set_ylim(0, 115)
    ax1.set_xticks(xs, labels, fontsize=9)
    ax1.set_ylabel('Rb + Pantoea UFBoot support (%)')
    ax1.set_title('A   Topological robustness', loc='left', fontsize=12, fontweight='bold')
    for x, v in zip(xs, boots):
        ax1.text(x, v + 2, f'{v:g}%', ha='center', fontsize=10, fontweight='bold') if v is not None else None
    ax2.bar(xs, logp, color=['#8094aa', '#167b78', '#167b78'], width=0.58)
    ax2.axhline(-math.log10(0.05), color='#b45449', lw=1.2, ls='--')
    ax2.text(2.48, -math.log10(0.05) + 0.5, 'p = 0.05', color='#b45449', fontsize=8, ha='right')
    ax2.set_xticks(xs, labels, fontsize=9)
    ax2.set_ylabel('−log₁₀(AU p-value for Rb + Rc)')
    ax2.set_title('B   Rejection of species-compatible topology', loc='left', fontsize=12, fontweight='bold')
    for x, v, p in zip(xs, logp, pvals):
        ax2.text(x, v + 0.7, f'p = {p:.2g}', ha='center', fontsize=9)
    ax2.set_ylim(0, max(logp) * 1.17)
    for ax in (ax1, ax2):
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#e7eaee', lw=0.6)
        ax.set_axisbelow(True)
    fig.text(0.5, 0.015,
             'AU tests compare candidate phylogenetic topologies (10,000 RELL resamplings); p-values are not probabilities of HGT.',
             ha='center', fontsize=9, color='#4b5563')
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.23, top=0.88)
    save_figure(fig, out, 'Figure3_phylogenetic_robustness')


def supplementary_tree_triptych(tree_map, out):
    fig, axes = plt.subplots(1, 3, figsize=(29, 10), gridspec_kw={'wspace': 0.22})
    for ax, (label, tree) in zip(axes, tree_map.items()):
        s = split_support(tree, focal_set(tree, 'RB', 'PA'))
        draw_cladogram(ax, tree, label, ('RB', 'PA'), s, font_size=6.5)
    fig.text(0.5, 0.025, '23 taxa per tree | Cladograms, not to scale | Rooted on D. solani for display only.',
             ha='center', fontsize=9, color='#4b5563')
    fig.subplots_adjust(left=0.02, right=0.985, bottom=0.08, top=0.89)
    save_figure(fig, out, 'FigureS1_three_sol_tree_robustness')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--root', type=Path, default=Path('/scratch/al98750/Roux/07_sol_HGT/07_phylogeny'))
    parser.add_argument('--out', type=Path, default=Path('/scratch/al98750/Roux/07_sol_HGT/10_manuscript_figures'))
    parser.add_argument('--summary-only', action='store_true', help='Produce numerical figures and tables without tree files (for previews).')
    args = parser.parse_args()
    root = args.root
    sol11 = root / 'sol11'
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 9, 'pdf.fonttype': 42,
        'svg.fonttype': 'none', 'axes.labelcolor': '#2b3139', 'text.color': '#2b3139',
        'axes.edgecolor': '#a5aeb8', 'figure.facecolor': 'white',
    })
    print('Reading completed results (no phylogenetic jobs will run)...')
    gene_rows = read_tsv(sol11 / 'robustness/gene_tree_focal_splits.tsv')
    robustness = read_tsv(sol11 / 'robustness/robustness_summary.tsv')
    if len(gene_rows) != 11:
        raise ValueError(f'Expected 11 genes, got {len(gene_rows)}')
    if len(robustness) != 2:
        raise ValueError(f'Expected two robustness analyses, got {len(robustness)}')
    robust_by_name = {r['analysis']: r for r in robustness}
    for expected in ('partitioned_11', 'partitioned_no_solA_10'):
        if expected not in robust_by_name:
            raise ValueError(f'Missing robustness analysis: {expected}')
    original_au = au_results(sol11 / 'targeted_AU/focal_AU_10k.iqtree')
    original_row = {
        'analysis': 'unpartitioned_11', 'taxa': '23', 'aligned_aa': '12987',
        'RB_PA_bootstrap': 'NA', 'RB_RC_bootstrap': 'NA',
        'AU_unconstrained': f"{original_au['unconstrained']:.10g}",
        'AU_RB_RC': f"{original_au['RB_RC']:.10g}",
        'AU_RB_PA': f"{original_au['RB_PA']:.10g}",
    }
    rows = [original_row, robust_by_name['partitioned_11'], robust_by_name['partitioned_no_solA_10']]
    if not args.summary_only:
        tree_files = {
            'species_core': root / 'species_core.treefile',
            'unpartitioned_11': sol11 / 'sol11_ML.treefile',
            'partitioned_11': sol11 / 'robustness/partitioned11/sol11_partitioned_ML.treefile',
            'no_solA_10': sol11 / 'robustness/noSolA_partitioned10/sol10_partitioned_ML.treefile',
        }
        trees = {name: Phylo.read(path, 'newick') for name, path in tree_files.items()}
        print('Validating taxon labels and group membership...')
        validate_taxa(trees)
        original_row['RB_PA_bootstrap'] = str(split_support(trees['unpartitioned_11'], focal_set(trees['unpartitioned_11'], 'RB', 'PA')) or 'absent')
        original_row['RB_RC_bootstrap'] = str(split_support(trees['unpartitioned_11'], focal_set(trees['unpartitioned_11'], 'RB', 'RC')) or 'absent')
        print('Creating Figure 1 (species vs partitioned sol trees)...')
        species_support, sol_support = figure1(trees['species_core'], trees['partitioned_11'], out)
        print(f'  Core species Rb+Rc split UFBoot: {species_support}')
        print(f'  Partitioned sol Rb+Pantoea split UFBoot: {sol_support}')
        supplementary_tree_triptych({
            'A  Unpartitioned 11 genes': trees['unpartitioned_11'],
            'B  Partitioned 11 genes': trees['partitioned_11'],
            'C  Partitioned 10 genes (no solA)': trees['no_solA_10'],
        }, out)
    else:
        # Preview mode uses numerical summaries; no tree bootstrap value can be inferred.
        original_row['RB_PA_bootstrap'] = 'NA'
    print('Creating Figure 2 (individual gene support)...')
    figure2(gene_rows, out)
    print('Creating Figure 3 (robustness comparison)...')
    if args.summary_only:
        # No tree file: use the reported original unpartitioned bootstrap only if verified externally.
        # Here plot two analyses, with the first marked as unavailable rather than fabricated.
        print('  NOTE: Original 11-gene bootstrap is unavailable in summary-only mode; marked unverified (not 0).')
    figure3(rows, out)
    table1_cols = list(gene_rows[0].keys())
    table2_cols = ['analysis', 'taxa', 'aligned_aa', 'RB_PA_bootstrap', 'RB_RC_bootstrap',
                   'AU_unconstrained', 'AU_RB_RC', 'AU_RB_PA']
    for sep, ext in [('\t', 'tsv'), (',', 'csv')]:
        write_table(out / f'Table1_individual_gene_trees.{ext}', gene_rows, table1_cols, sep)
        write_table(out / f'Table2_targeted_AU_robustness.{ext}', rows, table2_cols, sep)
    readme = (
        'SOL HGT MANUSCRIPT FIGURES\n\n'
        'Figure1_species_vs_partitioned_sol_tree: Species-core versus partitioned 11-gene sol tree.\n'
        '  Visual root placed on Dickeya solani, not inferred to be a biological outgroup.\n'
        '  Cladogram branches are NOT proportional to evolutionary distances.\n'
        'Figure2_individual_gene_support: UFBoot support for the RB+Pantoea split and alignment lengths.\n'
        'Figure3_phylogenetic_robustness: UFBoot and -log10(AU p-value).\n'
        'FigureS1_three_sol_tree_robustness: Full 23-tip trees for the three sol analyses.\n'
        'Table1: Directly transcribed from gene_tree_focal_splits.tsv.\n'
        'Table2: Original AU values parsed from IQ-TREE report and partitioned values from robustness_summary.tsv.\n'
        'NOTE: Absence of a recovered split does not mean bootstrap = 0.\n'
        'NOTE: The linked sol genes are not independent transfer events.\n'
        'NOTE: Verify inferred topology constraints and sample identity before submission.\n'
    )
    (out / 'README_figures.txt').write_text(readme)
    print(f'\nFinished. Output: {out}')
    for f in sorted(out.iterdir()):
        print(' ', f.name)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(1)
