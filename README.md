# Gyne-DEG Atlas

Tissue-stratified differential-expression meta-analysis across **PCOS, endometriosis and
adenomyosis**, built from uniformly pre-quantified public datasets. Designed to run on
**8 GB RAM** machines (no genome alignment; pre-computed counts only).

## What this repo does
1. Pulls analysis-ready counts (recount3 / ARCHS4) or GEO processed matrices.
2. Runs per-dataset DEG (DESeq2 for RNA-seq, limma for microarray).
3. Meta-analyses DEGs **within each (condition x tissue) stratum** (random-effects + rank aggregation).
4. Functional enrichment + PPI hub genes.
5. A cross-dataset-validated ML gene panel.
6. Feeds robust hub targets into network pharmacology / docking (separate repos).

## Golden rules (do not break)
- **Never pool across tissues.** One analysis per (condition x tissue) stratum.
- **One primary tissue per dataset.** Auto tissue tags overlap; fix in the sample sheet.
- **RNA-seq and microarray are analysed separately**, then combined at the results level.
- **Human only.** Drop non-human / mixed-species series.
- See `docs/tissue_rules.md` before adding any dataset.


## Script status (implemented)
All pipeline scripts are implemented and ready to run in the conda envs:
- `fetch_expr.R` - recount3 (RNA-seq) / GEO series matrix (microarray) / local Salmon counts. Writes
  `<id>.expr.rds` plus `<id>.expr.csv.gz` and `<id>.meta.csv` (the latter two feed the ML step).
- `deg_deseq2.R` / `deg_limma.R` - per-dataset DEG in the standardized schema.
- `meta_combine.R` - RRA (default) or random-effects (metafor); falls back to a flagged DESCRIPTIVE
  result when a stratum has fewer than `meta.min_datasets`.
- `enrich_hubs.R` - clusterProfiler GO:BP + KEGG, and STRING v12 PPI hub genes.
- `ml_panel.py` - leave-one-DATASET-out elastic-net panel; RNA-seq and microarray kept separate;
  per-dataset z-scoring removes batch scale. (Tested on synthetic data.)

Needs internet at runtime: recount3/GEO downloads (fetch), KEGG (enrich), and the first STRING call
(downloads the human DB once). None of these run on the 8 GB constraint as alignment - they are small
metadata/DB fetches.


> **Windows users:** do not use conda for the R steps - it fails to solve Bioconductor.
> Install R from CRAN, run `Rscript install_R_packages.R` once, and follow
> `docs/windows_no_conda.md`. Run Snakemake WITHOUT `--use-conda`.

## Quick start
```bash
mamba env create -f environment.yml          # base (snakemake)
conda activate gyne-deg
# per-step envs are created automatically by snakemake --use-conda
snakemake -n                                  # dry run (shows plan)
snakemake --use-conda --cores 2 -p            # 8 GB: keep cores low
```

## Layout
```
config/      sample sheet, strata table, parameters
workflow/    Snakefile + rules/ + scripts/
docs/        DEG schema, tissue rules, contributing
data/        inputs (gitignored)
results/     outputs (gitignored)
```

## 8 GB survival notes
- Use 1-2 cores; counts matrices are small, but R can spike. Call `gc()` often.
- For WGCNA, filter to the top ~5000 variable genes (see `scripts/wgcna.R` TODO).
- Do NOT download FASTQ on these machines. Alignment (if ever needed) runs only on the 64 GB workstation.
