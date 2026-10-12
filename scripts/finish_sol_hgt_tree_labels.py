#!/usr/bin/env python3
"""Complete verified NCBI strain labels and create a matched 23-taxon species tree.

Requires Python 3.9+ and Biopython. NCBI metadata requires HTTPS access; if
unavailable, use --skip-ncbi to generate the matched tree from current labels.
Does not change IQ-TREE originals. Never treats an isolate ID as a verified strain.
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen
from Bio import Phylo

ROOT_DEFAULT = Path('/scratch/al98750/Roux/07_sol_HGT')
BAD_VALUES = {'', 'na', 'n/a', 'none', 'null', 'unknown', 'not available', 'not applicable', 'missing', 'unspecified'}


def read_tsv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, columns):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def nested(obj, *keys):
    for key in keys:
        if not isinstance(obj, dict):
            return ''
        obj = obj.get(key, '')
    return obj if isinstance(obj, str) else ''


def clean(value):
    v = ' '.join((value or '').split())
    return '' if v.casefold() in BAD_VALUES else v


def acc_base(a):
    return a.split('.')[0]


def query_ncbi(accessions):
    """Fetch reports in small batches; retry failed batches individually."""
    urlbase = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'
    headers = {'Accept': 'application/json', 'User-Agent': 'sol-hgt-strain-metadata/1.0'}
    if os.getenv('NCBI_API_KEY'):
        headers['api-key'] = os.environ['NCBI_API_KEY']
    found = {}
    errors = []

    def get_group(group):
        path = quote(','.join(group), safe=',._')
        req = Request(urlbase + path + '/dataset_report?page_size=100', headers=headers)
        with urlopen(req, timeout=45) as res:
            payload = json.load(res)
        reports = payload.get('reports', [])
        if not isinstance(reports, list):
            raise ValueError('NCBI response did not contain a reports list')
        return reports

    for i in range(0, len(accessions), 8):
        group = accessions[i:i + 8]
        try:
            reports = get_group(group)
        except (HTTPError, URLError, TimeoutError, ValueError, json.JSONDecodeError) as e:
            print(f'WARNING: NCBI batch lookup failed ({type(e).__name__}: {e}); trying one by one', file=sys.stderr)
            reports = []
            for a in group:
                try:
                    reports.extend(get_group([a]))
                except (HTTPError, URLError, TimeoutError, ValueError, json.JSONDecodeError) as e2:
                    errors.append(f'{a}: {type(e2).__name__}: {e2}')
                time.sleep(0.25)
        for r in reports:
            a = r.get('accession', '')
            if a:
                found[a] = r
        print(f'NCBI batches: {min(i + 8, len(accessions))}/{len(accessions)} accessions attempted')
        time.sleep(0.3)
    return found, errors


def get_report(accession, found):
    if accession in found:
        return found[accession], 'exact'
    # Unversioned input accessions may resolve to a single versioned accession.
    matches = [r for a, r in found.items() if acc_base(a) == acc_base(accession)]
    if len(matches) == 1:
        return matches[0], 'unversioned_match'
    return None, 'not_found' if not matches else 'ambiguous_version'


def extract_record(acc, manifest_organism, report, match_type):
    if not report:
        return dict(accession=acc, manifest_organism=manifest_organism,
                    ncbi_accession='', ncbi_organism='', ncbi_strain='', ncbi_isolate='',
                    biosample='', status=match_type)
    org = report.get('organism') or {}
    info = report.get('assembly_info') or {}
    bio = info.get('biosample') or {}
    if not isinstance(bio, dict):
        bio = {}
    name = clean(nested(org, 'organism_name'))
    strain = clean(nested(org, 'infraspecific_names', 'strain')) or clean(nested(bio, 'strain')) or clean(nested(bio, 'description', 'organism', 'infraspecific_names', 'strain'))
    isolate = clean(nested(org, 'infraspecific_names', 'isolate')) or clean(nested(bio, 'isolate')) or clean(nested(bio, 'description', 'organism', 'infraspecific_names', 'isolate'))
    # Only auto-fill if NCBI taxonomic name is compatible with the manifest.
    prefix = ' '.join(manifest_organism.split()[:2]).casefold()
    taxon_ok = not name or name.casefold().startswith(prefix + ' ') or name.casefold() == prefix
    status = ('matched' if taxon_ok else 'TAXON_MISMATCH_REVIEW')
    return dict(accession=acc, manifest_organism=manifest_organism,
                ncbi_accession=report.get('accession', ''), ncbi_organism=name,
                ncbi_strain=strain, ncbi_isolate=isolate,
                biosample=clean(nested(info, 'biosample', 'accession')),
                status=status)


def prune_matched_tree(out):
    core = out / '01_species_core_32taxa_labeled.nwk'
    sol = out / '02_sol11_partitioned_23taxa_labeled.nwk'
    coretree = Phylo.read(core, 'newick')
    soltree = Phylo.read(sol, 'newick')
    sol_tips = {x.name for x in soltree.get_terminals()}
    core_tips = {x.name for x in coretree.get_terminals()}
    if len(sol_tips) != 23 or len(core_tips) != 32 or not sol_tips.issubset(core_tips):
        raise SystemExit(f'ERROR: taxa mismatch: species={len(core_tips)}, sol={len(sol_tips)}, missing={sol_tips - core_tips}')
    for tip in list(coretree.get_terminals()):
        if tip.name not in sol_tips:
            coretree.prune(tip)
    target = out / '05_species_core_matched23_labeled.nwk'
    Phylo.write(coretree, target, 'newick')
    final = Phylo.read(target, 'newick')
    final_tips = [x.name for x in final.get_terminals()]
    if len(final_tips) != 23 or set(final_tips) != sol_tips:
        raise SystemExit('ERROR: pruned tree QC failed')
    write_tsv(out / 'matched23_tree_QC.tsv', [
        {'file': core.name, 'taxa': 32, 'purpose': 'full original-inference species tree'},
        {'file': sol.name, 'taxa': 23, 'purpose': 'main sol tree'},
        {'file': target.name, 'taxa': 23, 'purpose': 'display-only pruned species tree; support from 32-taxon analysis'},
    ], ['file', 'taxa', 'purpose'])
    print(f'PASS: matched species tree has 23 taxa: {target}')
    print('NOTE: pruning is for figure display only; internal support was not recalculated.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT_DEFAULT)
    p.add_argument('--label-script', type=Path, default=Path.home() / 'Rb_Solanimycin/scripts/label_sol_hgt_trees.py')
    p.add_argument('--skip-ncbi', action='store_true', help='Skip online metadata; just relabel and make matched 23-taxon species tree')
    args = p.parse_args()
    root = args.root
    out = root / 'Geneious_trees_labeled'
    out.mkdir(parents=True, exist_ok=True)
    manifest = read_tsv(root / '05_results/genome_manifest.tsv')
    accession_to_org = {r['accession'].strip(): r['organism'].strip() for r in manifest}
    if len(accession_to_org) != 32:
        raise SystemExit(f'ERROR: expected 32 distinct manifest accessions, got {len(accession_to_org)}')
    overrides_path = out / 'strain_overrides.tsv'
    if overrides_path.exists():
        existing = {r['accession'].strip(): (r.get('strain') or '').strip() for r in read_tsv(overrides_path)}
    else:
        existing = {}
    additions = 0
    if not args.skip_ncbi:
        found, errors = query_ncbi(list(accession_to_org))
        records = []
        for a, org in accession_to_org.items():
            report, match_type = get_report(a, found)
            r = extract_record(a, org, report, match_type)
            strain = r['ncbi_strain']
            # Preserve user-entered strain labels. Never promote an isolate to a strain.
            if r['status'] == 'matched' and strain and not existing.get(a):
                if strain.casefold() not in org.casefold():
                    existing[a] = strain
                    additions += 1
            r['strain_in_label'] = existing.get(a, '')
            records.append(r)
        write_tsv(out / 'ncbi_assembly_metadata.tsv', records,
                  ['accession', 'manifest_organism', 'ncbi_accession', 'ncbi_organism',
                   'ncbi_strain', 'ncbi_isolate', 'biosample', 'status', 'strain_in_label'])
        if errors:
            (out / 'ncbi_lookup_errors.txt').write_text('\n'.join(errors) + '\n')
            print(f'WARNING: {len(errors)} individual NCBI requests failed; see ncbi_lookup_errors.txt')
        print(f'NCBI metadata matched: {sum(r["status"] == "matched" for r in records)}/32')
        print(f'New verified strain labels: {additions}')
        print('Metadata for review: ' + str(out / 'ncbi_assembly_metadata.tsv'))
    write_tsv(overrides_path,
              [{'accession': a, 'strain': existing.get(a, '')} for a in sorted(accession_to_org)],
              ['accession', 'strain'])
    if not args.label_script.is_file():
        raise SystemExit(f'ERROR: existing label script not found: {args.label_script}')
    subprocess.run([sys.executable, str(args.label_script), '--root', str(root)], check=True)
    prune_matched_tree(out)
    print('\nNEXT: review ncbi_assembly_metadata.tsv for missing strain names, then import 05 and 02 into Geneious.')


if __name__ == '__main__':
    main()
