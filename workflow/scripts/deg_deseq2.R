#!/usr/bin/env Rscript
# Per-dataset DEG via DESeq2. Writes the standardized schema (see docs/DEG_output_schema.md).
suppressMessages({library(optparse); library(DESeq2); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--expr"), make_option("--group"), make_option("--case"),
  make_option("--control"), make_option("--padj",type="double",default=0.05),
  make_option("--lfc",type="double",default=1.0), make_option("--out"))))
x <- readRDS(opt$expr); counts <- round(as.matrix(x$counts)); meta <- x$meta
meta[[opt$group]] <- factor(meta[[opt$group]], levels=c(opt$control, opt$case))
keep <- rowSums(counts >= 10) >= max(2, floor(0.2*ncol(counts)))   # low-count filter
dds <- DESeqDataSetFromMatrix(counts[keep,], meta, as.formula(paste("~", opt$group)))
dds <- DESeq(dds); res <- as.data.frame(results(dds))
res$gene <- rownames(res)
out <- data.table(gene=res$gene, log2FC=res$log2FoldChange, p=res$pvalue,
                  padj=res$padj, n=ncol(counts))
dir.create(dirname(opt$out), showWarnings=FALSE, recursive=TRUE)
fwrite(out, opt$out, sep="\t")
gc()
