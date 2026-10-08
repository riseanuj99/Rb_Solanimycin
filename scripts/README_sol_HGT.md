# 20GA0316 *sol* evolutionary analysis — Sapelo2

## Goal

Run six analyses: (1) core-proteome species tree; (2) *sol* phylogeny; (3) species-versus-*sol* tree comparison and AU test; (4) *sol* distribution across *Rouxiella*; (5) cluster synteny and flanks; (6) candidate insertion sites, mobile-element annotations and GC composition.

**Important:** Matilla et al. (2022) already demonstrated evidence for intergeneric horizontal transfer of *sol*. This pipeline tests the specific history and genomic integration of the 20GA0316 locus. It does not assume an HGT result.

## Install once on Sapelo2

```bash
mkdir -p /scratch/al98750/Roux/logs /scratch/al98750/Roux/envs
source "$(conda info --base)/etc/profile.d/conda.sh"
mamba create -y -p /scratch/al98750/Roux/envs/sol_hgt -c conda-forge -c bioconda \
  python=3.11 biopython ncbi-datasets-cli prodigal blast mafft trimal iqtree orthofinder diamond fasttree
```

If `mamba` is unavailable, substitute `conda create` (same arguments). Installing `clinker` is optional for visualization; it is not needed for data discovery or tree inference. No packages are installed into your existing `roux` environment.

## Put scripts into your laptop GitHub checkout

Copy all three scripts (`03_sol_HGT.py`, `03_sol_HGT_discovery.sh`, `03_sol_HGT_phylogeny.sh`) to `Rb_Solanimycin/scripts/` on your laptop, then commit and push:

```bash
git add scripts/03_sol_HGT.py scripts/03_sol_HGT_discovery.sh scripts/03_sol_HGT_phylogeny.sh
git commit -m "Add sol cluster evolutionary analysis"
git push
```

## Submit on Sapelo2

```bash
cd ~/Rb_Solanimycin
git pull
mkdir -p /scratch/al98750/Roux/logs
sbatch scripts/03_sol_HGT_discovery.sh
```

Wait until the discovery job finishes successfully and review `07_sol_HGT/05_results/sol_loci.tsv`, especially the positive controls, before submitting phylogeny:

```bash
sbatch scripts/03_sol_HGT_phylogeny.sh
```

## Key outputs

All under `/scratch/al98750/Roux/07_sol_HGT/`:

- `05_results/genome_manifest.tsv`: 17 local *R. badensis* assemblies plus up to three additional assemblies per other *Rouxiella* species and published *sol*-positive genera.
- `02_reference/sol_reference.json`: validated original locus tags `LLR01_11590` through `LLR01_11530` and exact contig boundaries.
- `05_results/sol_loci.tsv`, `sol_gene_matrix.tsv`, `sol_gene_hits.tsv`: candidate *sol* presence/absence, hits and coordinates.
- `06_neighborhoods/*_sol_25kb.gb`: 25-kb flanks around detected *sol* clusters, with CDS translations for Clinker.
- `05_results/flank_anchor_comparison.tsv`, `integration_site_candidates.tsv`: potential conserved/empty integration-site comparisons (preliminary candidates, not definitive boundaries).
- `05_results/mobile_element_candidates.tsv`: annotated mobile-element/tRNA candidates around detected clusters; missing annotation is marked explicitly.
- `08_composition/GC_permutation.json`: same-length genomic-window GC comparison and exploratory empirical p value.
- `07_phylogeny/species_core.treefile`: maximum-likelihood species tree from up to 200 single-copy core protein families.
- `07_phylogeny/sol_cluster.treefile`: ML tree from a curated set of *sol* tailoring/transport genes, excluding modular PKS/NRPS genes and lineage-specific *solM*.
- `07_phylogeny/species_pruned_to_sol_taxa.nwk`, `tree_comparison.json`, `topology_AU_test.iqtree`: same-taxon tree comparison and statistical topology test.

## Optional Clinker visualization

After installing `clinker` in a suitable environment, check `clinker --help`, then plot selected neighborhood GenBank files; e.g.:

```bash
cd /scratch/al98750/Roux/07_sol_HGT/06_neighborhoods
clinker GCF_020740305.1_sol_25kb.gb GCA_000365285.1_sol_25kb.gb -p sol_comparison.html
```

The *sol* calls are initial homology/synteny candidates. For a manuscript, inspect false-positive PKS/NRPS matches, check reference gene orthology and domain architecture, and verify putative insertion junctions with nucleotide-level alignment and contig-boundary checks. No detected homolog is **not** automatically equivalent to confirmed biological absence. A phylogenetic discordance is not by itself proof of HGT.

## Checkpoints

`07_sol_HGT/.state/*.done` plus per-genome Prodigal/BLAST markers. Successful stages skip on resubmission. Do not delete checkpoint files unless you are deliberately invalidating upstream inputs. If you change reference queries or genome selection, invalidate dependent stages before rerunning.
