
#!/usr/bin/env Rscript
# One-shot installer for the WHOLE pipeline using a normal CRAN R install.
# Use this on Windows (and it works on Mac/Linux too) to AVOID the conda/bioconda
# solver problems. Bioconductor ships precompiled Windows binaries, so this just
# downloads them - no compiling, no environment solving.
#
# Run once per machine:   Rscript install_R_packages.R
options(repos = c(CRAN = "https://cloud.r-project.org"))
cran <- c("optparse","data.table","metafor","RobustRankAggreg","igraph","BiocManager")
for (p in cran) if (!requireNamespace(p, quietly=TRUE)) {
  message("Installing CRAN: ", p); install.packages(p)
}
bioc <- c("SummarizedExperiment","Biobase","DESeq2","edgeR","limma","GEOquery",
          "recount3","clusterProfiler","org.Hs.eg.db","STRINGdb")
message("Installing Bioconductor packages (precompiled binaries) ...")
BiocManager::install(bioc, update = FALSE, ask = FALSE)

need <- c(cran, bioc)
ok <- vapply(need, requireNamespace, logical(1), quietly = TRUE)
cat("\n==== install check ====\n")
for (i in seq_along(need)) cat(sprintf("  %-22s %s\n", need[i], if (ok[i]) "OK" else "MISSING"))
if (all(ok)) cat("\nALL GOOD - you can run the R scripts now (see docs/windows_no_conda.md).\n") else
  cat("\nSTILL MISSING:", paste(need[!ok], collapse=", "), "\n  Re-run, or install those individually.\n")
