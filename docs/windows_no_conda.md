
# Running on Windows without conda (recommended for the 8 GB team)

conda + Bioconductor often fails to solve on Windows. Skip it. Install R directly -
Bioconductor provides precompiled Windows binaries, so there is no solver and no
compiling. This is the reliable path.

## One-time setup (per machine)
1. Install R from CRAN:  https://cran.r-project.org/bin/windows/base/
   (Optional but recommended: install Rtools from the same site.)
2. Install all pipeline packages in one shot (run from the repo folder):
   ```
   Rscript install_R_packages.R
   ```
   This downloads DESeq2, limma, GEOquery, recount3, clusterProfiler, STRINGdb, etc.
   It prints OK / MISSING for each. Re-run if anything is MISSING.

## Smoke test WITHOUT conda or Snakemake (fastest way to get unblocked)
Run the two R scripts directly. Replace the labels in step 3 with the exact text you
read from the metadata file in step 2.

1. Fetch the expression + metadata for GSE5090:
   ```
   Rscript workflow/scripts/fetch_expr.R --dataset GSE5090 --technique Microarray --source GEO_matrix --out data/expr/GSE5090.expr.rds
   ```
2. Open data/expr/GSE5090.meta.csv. Find the column that splits PCOS vs control and the
   EXACT label text. (fetch_expr.R also prints the column names in the console.)
3. Build the DEG table (quote the values exactly as they appear):
   ```
   Rscript workflow/scripts/deg_limma.R --expr data/expr/GSE5090.expr.rds --group "disease state:ch1" --case "PCOS" --control "control" --out results/deg/GSE5090.deg.microarray.tsv
   ```
4. Open results/deg/GSE5090.deg.microarray.tsv. Confirm columns gene, log2FC, p, padj, n,
   real gene symbols, and SOME (not all, not none) genes with padj < 0.05.

If it errors, copy the FULL console output for a fix.

## Running the whole pipeline without conda
Once packages are installed, run Snakemake WITHOUT --use-conda (it will use your system R):
```
snakemake --cores 1 -p results/deg/GSE5090.deg.microarray.tsv
```
Do NOT pass --use-conda on Windows. (--use-conda is only for Mac/Linux/Colab, where the
env YAMLs in envs/ solve fine, including the new envs/fetch.yaml.)

## RNA-seq path
Repeat the smoke test once with a small RNA-seq dataset (source = recount3), building
results/deg/<GSE>.deg.tsv and using deg_deseq2.R. recount3 downloads precomputed counts;
no alignment, so it is fine on 8 GB.
