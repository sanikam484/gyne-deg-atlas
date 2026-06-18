#!/usr/bin/env Rscript
# Combine per-dataset DEGs WITHIN one (condition x tissue) stratum.
#   method = RRA  -> Robust Rank Aggregation (up & down separately)
#   method = REM  -> random-effects effect-size meta-analysis (metafor), SE derived from p
# If the stratum has < min_datasets clean datasets, writes a DESCRIPTIVE result (clearly flagged),
# never a meta-analysis.
suppressMessages({library(optparse); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--stratum"), make_option("--method", default="RRA"),
  make_option("--min_datasets", type="integer", default=4L),
  make_option("--deg_dir", default="results/deg"),
  make_option("--samplesheet"), make_option("--out"))))

norm_id <- function(cond, tis) gsub("/", "-", gsub(" ", "_", paste0(cond, "__", tis)))
ss <- fread(opt$samplesheet)
ss[, sid := norm_id(condition, tissue_primary)]
if (!"notes" %in% names(ss)) ss[, notes := ""]
ids <- ss[sid == opt$stratum & (is.na(notes) | !grepl("DROP", notes, ignore.case=TRUE)), dataset_id]
ids <- unique(ids)
if (length(ids) == 0) { fwrite(data.table(note="no datasets in stratum"), opt$out, sep="\t"); quit(save="no") }

read_deg <- function(d){
  for (f in c(file.path(opt$deg_dir, paste0(d, ".deg.tsv")),
              file.path(opt$deg_dir, paste0(d, ".deg.microarray.tsv")))) {
    if (file.exists(f)) { x <- fread(f); x[, dataset := d]; return(x) }
  }
  message("[meta] missing DEG for ", d); NULL
}
degs <- rbindlist(lapply(ids, read_deg), use.names=TRUE, fill=TRUE)
degs <- degs[!is.na(gene) & gene != ""]
nds  <- length(unique(degs$dataset))
message(sprintf("[meta] stratum %s: %d datasets with DEGs.", opt$stratum, nds))

write_out <- function(dt){ fwrite(dt, opt$out, sep="\t"); message("[meta] wrote ", opt$out) }

## ---- descriptive fallback (too few datasets) ----
if (nds < opt$min_datasets) {
  sig <- degs[!is.na(padj) & padj < 0.05]
  if (nrow(sig) == 0) sig <- degs
  desc <- sig[, .(n_datasets = uniqueN(dataset),
                  mean_log2FC = mean(log2FC, na.rm=TRUE),
                  min_padj    = min(padj, na.rm=TRUE)), by = gene]
  desc[, direction := ifelse(mean_log2FC > 0, "up", "down")]
  desc[, note := sprintf("DESCRIPTIVE only (%d dataset(s) < min_datasets=%d); NOT a meta-analysis",
                         nds, opt$min_datasets)]
  setorder(desc, min_padj)
  write_out(desc); quit(save="no")
}

## ---- Robust Rank Aggregation ----
if (toupper(opt$method) == "RRA") {
  suppressMessages(library(RobustRankAggreg))
  sig <- degs[!is.na(padj) & padj < 0.05]
  agg_dir <- function(dir){
    parts <- split(sig, sig$dataset)
    lists <- lapply(parts, function(d){
      d <- if (dir == "up") d[log2FC > 0] else d[log2FC < 0]
      d <- d[order(-abs(log2FC))]; as.character(d$gene)
    })
    lists <- lists[vapply(lists, length, 1L) > 0]
    if (length(lists) < 2) return(data.table(gene=character(), meta_p=numeric(), direction=character()))
    r <- as.data.table(aggregateRanks(glist = lists, full = TRUE))
    setnames(r, c("Name", "Score"), c("gene", "meta_p")); r[, direction := dir]; r[]
  }
  res <- rbind(agg_dir("up"), agg_dir("down"))
  if (nrow(res) == 0) { write_out(data.table(note="no significant genes for RRA")); quit(save="no") }
  res[, meta_padj := p.adjust(meta_p, "BH")][, n_datasets := nds]
  setorder(res, meta_p); res <- res[!duplicated(gene)]
  write_out(res[, .(gene, meta_p, meta_padj, direction, n_datasets)]); quit(save="no")
}

## ---- random-effects (metafor) ----
if (toupper(opt$method) == "REM") {
  suppressMessages(library(metafor))
  d <- degs[!is.na(log2FC) & !is.na(p)]
  d[, z := qnorm(pmin(pmax(p, 1e-300), 1) / 2, lower.tail = FALSE)]   # |z| from two-sided p
  d[, se := fifelse(z > 0, abs(log2FC) / z, NA_real_)]
  d <- d[is.finite(se) & se > 0]
  genes <- d[, .N, by = gene][N >= 2, gene]
  out <- rbindlist(lapply(genes, function(g){
    s <- d[gene == g]
    fit <- tryCatch(rma(yi = s$log2FC, sei = s$se, method = "REML"),
             error = function(e) tryCatch(rma(yi = s$log2FC, sei = s$se, method = "DL"),
                                          error = function(e) NULL))
    if (is.null(fit)) return(NULL)
    data.table(gene=g, meta_log2FC=as.numeric(fit$b), meta_se=fit$se,
               meta_p=fit$pval, I2=fit$I2, n_datasets=nrow(s))
  }), fill=TRUE)
  if (nrow(out) == 0) { write_out(data.table(note="no genes shared across >=2 datasets")); quit(save="no") }
  out[, meta_padj := p.adjust(meta_p, "BH")][, direction := ifelse(meta_log2FC > 0, "up", "down")]
  setorder(out, meta_p)
  write_out(out); quit(save="no")
}
stop("Unknown --method (use RRA or REM).")
