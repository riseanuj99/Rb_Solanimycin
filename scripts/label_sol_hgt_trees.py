#!/usr/bin/env python3
"""Relabel existing IQ-TREE Newick tips using a verified genome manifest.

Input: /scratch/al98750/Roux/07_sol_HGT/05_results/genome_manifest.tsv
Optional: strain_overrides.tsv (accession<TAB>strain), created as a template.
Output: labeled .nwk trees for Geneious, plus mapping and QC tables.

Only tip tokens are replaced in the original Newick text. Branch lengths,
internal support labels, tree structure, and all other text remain unchanged.
"""
import argparse
import csv
import re
from pathlib import Path
from Bio import Phylo

ACCESSION = r"G(?:CA|CF)_\d{9}(?:\.\d+)?"
TIP_RE = re.compile(r"(?<=[(,])(\s*)(" + ACCESSION + r")(?=\s*[:),])")


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def label(organism, accession, strain):
    organism = organism.strip()
    strain = strain.strip()
    # Avoid duplicating strain designations already present in organism.
    if strain and strain.casefold() not in organism.casefold().split():
        organism += " " + strain
    return f"{organism} | {accession}"


def quote_newick(text):
    return "'" + text.replace("'", "''") + "'"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path("/scratch/al98750/Roux/07_sol_HGT"))
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    root = args.root
    out = args.out or root / "Geneious_trees_labeled"
    out.mkdir(parents=True, exist_ok=True)

    manifest_path = root / "05_results/genome_manifest.tsv"
    if not manifest_path.is_file():
        raise SystemExit(f"ERROR: manifest not found: {manifest_path}")
    rows = read_tsv(manifest_path)
    if not rows or not {"accession", "organism"}.issubset(rows[0]):
        raise SystemExit("ERROR: manifest needs accession and organism columns")
    metadata = {}
    for row in rows:
        acc = row["accession"].strip()
        org = row["organism"].strip()
        if acc in metadata and metadata[acc] != org:
            raise SystemExit(f"ERROR: conflicting organism names for {acc}")
        metadata[acc] = org

    overrides_path = out / "strain_overrides.tsv"
    if not overrides_path.exists():
        write_tsv(overrides_path,
                  [{"accession": acc, "strain": ""} for acc in sorted(metadata)],
                  ["accession", "strain"])
        print(f"Created optional strain template: {overrides_path}")
    overrides = {}
    for row in read_tsv(overrides_path):
        if "accession" not in row or "strain" not in row:
            raise SystemExit("ERROR: strain_overrides.tsv needs accession and strain columns")
        acc = row["accession"].strip()
        if acc not in metadata:
            raise SystemExit(f"ERROR: strain override accession not in manifest: {acc}")
        overrides[acc] = (row["strain"] or "").strip()

    phy = root / "07_phylogeny"
    files = [
        (phy / "species_core.treefile", "01_species_core_32taxa_labeled.nwk"),
        (phy / "sol11/robustness/partitioned11/sol11_partitioned_ML.treefile",
         "02_sol11_partitioned_23taxa_labeled.nwk"),
        (phy / "sol11/robustness/noSolA_partitioned10/sol10_partitioned_ML.treefile",
         "03_sol10_no_solA_23taxa_labeled.nwk"),
        (phy / "sol11/sol11_ML.treefile", "04_sol11_unpartitioned_23taxa_labeled.nwk"),
    ]
    files += [(p, f"gene_trees/{p.stem}_labeled.nwk")
              for p in sorted((phy / "sol11/alignments").glob("sol?_ML.treefile"))]

    used = set()
    results = []
    for source, relpath in files:
        if not source.is_file():
            raise SystemExit(f"ERROR: required tree missing: {source}")
        original = source.read_text()
        before = Phylo.read(str(source), "newick")
        tips = [t.name for t in before.get_terminals()]
        if len(tips) != len(set(tips)):
            raise SystemExit(f"ERROR: duplicate tip IDs in {source}")
        missing = set(tips) - set(metadata)
        if missing:
            raise SystemExit(f"ERROR: missing manifest metadata for {source}: {sorted(missing)}")
        counts = {}
        def replacement(match):
            whitespace, acc = match.groups()
            counts[acc] = counts.get(acc, 0) + 1
            used.add(acc)
            return whitespace + quote_newick(label(metadata[acc], acc, overrides.get(acc, "")))
        renamed = TIP_RE.sub(replacement, original)
        if set(counts) != set(tips) or any(count != 1 for count in counts.values()):
            raise SystemExit(f"ERROR: tip matching failed in {source}: "
                             f"found {len(counts)} of {len(tips)}; counts={counts}")
        target = out / relpath
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(renamed)
        after = Phylo.read(str(target), "newick")
        new_tips = [t.name for t in after.get_terminals()]
        expected = {label(metadata[acc], acc, overrides.get(acc, "")) for acc in tips}
        if set(new_tips) != expected or len(new_tips) != len(tips):
            target.unlink(missing_ok=True)
            raise SystemExit(f"ERROR: Geneious/Newick QC failed for {target}")
        results.append({"file": str(target), "taxa": len(tips), "status": "PASS"})
        print(f"PASS {len(tips):2d} taxa: {target}")

    write_tsv(out / "tip_label_mapping.tsv", [
        {"accession": acc, "organism": metadata[acc],
         "strain_override": overrides.get(acc, ""),
         "geneious_tip_label": label(metadata[acc], acc, overrides.get(acc, ""))}
        for acc in sorted(used)
    ], ["accession", "organism", "strain_override", "geneious_tip_label"])
    write_tsv(out / "tree_QC.tsv", results, ["file", "taxa", "status"])
    print(f"\nAll {len(results)} trees passed QC; {len(used)} unique accession labels.")
    print("No original tree files were modified.")
    print("Add verified strain names to strain_overrides.tsv and rerun to update labels.")

if __name__ == "__main__":
    main()
