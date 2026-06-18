#!/usr/bin/env Rscript
# Pull analysis-ready expression for ONE dataset and write:
#   <out>.expr.rds      list(counts = genes x samples, meta = samples x covariates)
#   <out sibling>.expr.csv.gz   genes x samples (for Python / inspection)
#   <out sibling>.meta.csv      samples x covariates
# Sources:
#   recount3          -> RNA-seq raw read counts (preferred; 8 GB friendly, one dataset at a time)
#   GEO_matrix        -> microarray processed series matrix (log2 intensities)
#   salmon_workstation-> read counts you produced on the 64 GB Mac (data/salmon/<GSE>_counts.csv)
#   archs4            -> not implemented in R here; export ARCHS4 to a counts CSV and use salmon_workstation
suppressMessages({library(optparse)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--dataset"), make_option("--technique", default=""),
  make_option("--source"),  make_option("--out"),
  make_option("--salmon_dir", default="data/salmon"))))

dataset <- opt$dataset; src <- tolower(opt$source)
out_rds  <- opt$out
out_expr <- sub("\\.expr\\.rds$", ".expr.csv.gz", out_rds)
out_meta <- sub("\\.expr\\.rds$", ".meta.csv",   out_rds)
dir.create(dirname(out_rds), showWarnings=FALSE, recursive=TRUE)

write_outputs <- function(counts, meta, integer_counts){
  counts <- as.matrix(counts)
  if (isTRUE(integer_counts)) storage.mode(counts) <- "integer"
  meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) || any(rownames(meta)=="")) rownames(meta) <- colnames(counts)
  saveRDS(list(counts=counts, meta=meta), out_rds)
  ec <- data.frame(gene=rownames(counts), counts, check.names=FALSE)
  gz <- gzfile(out_expr, "w"); write.csv(ec, gz, row.names=FALSE); close(gz)
  m  <- cbind(sample=rownames(meta), meta)
  write.csv(m, out_meta, row.names=FALSE)
  message(sprintf("[fetch] %s: %d genes x %d samples written.", dataset, nrow(counts), ncol(counts)))
  message("[fetch] meta columns: ", paste(colnames(meta), collapse=", "))
  message("[fetch] -> set group_column/case_label/control_label in the sample sheet to one of the meta columns above.")
}

if (src == "recount3") {
  suppressMessages({library(recount3); library(SummarizedExperiment)})
  ap  <- available_projects()
  sel <- ap[ap$project == dataset, , drop=FALSE]
  if (nrow(sel) == 0) stop(sprintf("recount3 has no project '%s'. Use GEO_matrix or salmon_workstation.", dataset))
  rse <- create_rse(sel[1, ])
  counts <- compute_read_counts(rse)                 # raw read counts for DESeq2
  rd  <- as.data.frame(rowData(rse))
  sym <- if ("gene_name" %in% colnames(rd)) rd$gene_name else rownames(counts)
  sym <- ifelse(is.na(sym) | sym == "", rownames(counts), sym)
  rownames(counts) <- make.unique(as.character(sym))
  write_outputs(round(counts), as.data.frame(colData(rse)), integer_counts=TRUE)

} else if (src == "geo_matrix") {
  suppressMessages({library(GEOquery); library(Biobase)})
  g <- getGEO(dataset, GSEMatrix=TRUE, getGPL=TRUE)
  if (length(g) == 0) stop(sprintf("GEO returned nothing for %s.", dataset))
  eset <- g[[1]]
  E <- exprs(eset)
  qx <- as.numeric(quantile(E, c(0, .25, .5, .75, .99, 1), na.rm=TRUE))
  if (qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0)) {   # looks linear -> log2
    E[E <= 0] <- NA; E <- log2(E)
  }
  fd <- as.data.frame(fData(eset))
  symcol <- intersect(c("Gene symbol","Gene Symbol","GENE_SYMBOL","Symbol","SYMBOL","gene_symbol"), colnames(fd))
  if (length(symcol) >= 1) {
    sym <- as.character(fd[[symcol[1]]]); sym <- sub(" ///.*", "", sym)
    keep <- !(is.na(sym) | sym == "" | sym == "---"); E <- E[keep, , drop=FALSE]; sym <- sym[keep]
    o <- order(rowMeans(E, na.rm=TRUE), decreasing=TRUE)   # keep highest-expressed probe per gene
    E <- E[o, , drop=FALSE]; sym <- sym[o]
    dup <- duplicated(sym); E <- E[!dup, , drop=FALSE]; rownames(E) <- sym[!dup]
  } else message("[fetch] No gene-symbol column in featureData; keeping probe IDs (map later).")
  write_outputs(E, as.data.frame(pData(eset)), integer_counts=FALSE)

} else if (src == "salmon_workstation") {
  cf <- file.path(opt$salmon_dir, paste0(dataset, "_counts.csv"))
  mf <- file.path(opt$salmon_dir, paste0(dataset, "_meta.csv"))
  if (!file.exists(cf)) stop(sprintf("Expected %s (produce it on the 64 GB workstation: Salmon + tximport -> gene counts CSV).", cf))
  counts <- as.matrix(read.csv(cf, row.names=1, check.names=FALSE))
  meta <- if (file.exists(mf)) read.csv(mf, row.names=1, check.names=FALSE) else data.frame(row.names=colnames(counts))
  write_outputs(round(counts), meta, integer_counts=TRUE)

} else if (src == "archs4") {
  stop("ARCHS4 not implemented in R here. Export the ARCHS4 counts to data/salmon/<GSE>_counts.csv and use data_source=salmon_workstation, or use recount3.")
} else {
  stop(sprintf("Unknown data_source '%s'. Use recount3 | GEO_matrix | salmon_workstation.", opt$source))
}
