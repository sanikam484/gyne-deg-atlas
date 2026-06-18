#!/usr/bin/env Rscript
# Functional enrichment (GO:BP + KEGG) and STRING PPI hub genes for one stratum's meta signature.
# First STRING use downloads a species database (needs internet); 9606 = human.
suppressMessages({library(optparse); library(data.table)})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--meta"), make_option("--out_enrich"), make_option("--out_hubs"),
  make_option("--padj", type="double", default=0.05),
  make_option("--top", type="integer", default=400L))))   # cap genes sent to STRING (8 GB friendly)

m <- fread(opt$meta)
if (!"gene" %in% names(m)) { fwrite(data.table(note="no gene column"), opt$out_enrich, sep="\t")
                             fwrite(data.table(note="no gene column"), opt$out_hubs, sep="\t"); quit(save="no") }
sigcol <- if ("meta_padj" %in% names(m)) "meta_padj" else if ("min_padj" %in% names(m)) "min_padj" else "meta_p"
sig <- m[is.finite(get(sigcol)) & get(sigcol) < opt$padj]
genes <- unique(as.character(sig$gene))
message(sprintf("[enrich] %d significant genes (%s < %g).", length(genes), sigcol, opt$padj))
if (length(genes) < 10) {
  fwrite(data.table(note=sprintf("only %d sig genes; enrichment skipped", length(genes))), opt$out_enrich, sep="\t")
  fwrite(data.table(note="too few genes for PPI"), opt$out_hubs, sep="\t"); quit(save="no")
}

## enrichment
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db)})
eg <- tryCatch(bitr(genes, "SYMBOL", "ENTREZID", OrgDb=org.Hs.eg.db), error=function(e) NULL)
enr <- data.table()
if (!is.null(eg) && nrow(eg) > 0) {
  ego <- tryCatch(enrichGO(eg$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP",
                  pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE),
                  error=function(e) NULL)
  ek  <- tryCatch(enrichKEGG(eg$ENTREZID, organism="hsa", pvalueCutoff=0.05), error=function(e) NULL)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) enr <- rbind(enr, data.table(source="GO:BP", as.data.table(as.data.frame(ego))), fill=TRUE)
  if (!is.null(ek)  && nrow(as.data.frame(ek))  > 0) enr <- rbind(enr, data.table(source="KEGG",  as.data.table(as.data.frame(ek))),  fill=TRUE)
}
if (nrow(enr) == 0) enr <- data.table(note="no enriched terms at p<0.05")
fwrite(enr, opt$out_enrich, sep="\t")

## STRING PPI hubs
hubs <- data.table(note="no PPI edges")
ok <- requireNamespace("STRINGdb", quietly=TRUE) && requireNamespace("igraph", quietly=TRUE)
if (ok) {
  suppressMessages({library(STRINGdb); library(igraph)})
  g_in <- head(genes, opt$top)
  sdb <- tryCatch(STRINGdb$new(version="12.0", species=9606, score_threshold=400, input_directory=tempdir()),
                  error=function(e) NULL)
  if (!is.null(sdb)) {
    mp <- sdb$map(data.frame(gene=g_in, stringsAsFactors=FALSE), "gene", removeUnmappedRows=TRUE)
    ints <- sdb$get_interactions(mp$STRING_id)
    if (!is.null(ints) && nrow(ints) > 0) {
      gr <- graph_from_data_frame(unique(ints[, c("from","to")]), directed=FALSE)
      deg <- igraph::degree(gr); btw <- igraph::betweenness(gr)
      id2sym <- setNames(mp$gene, mp$STRING_id)
      hubs <- data.table(STRING_id=names(deg), gene=unname(id2sym[names(deg)]),
                         degree=as.integer(deg), betweenness=as.numeric(btw[names(deg)]))
      setorder(hubs, -degree, -betweenness)
    }
  }
} else hubs <- data.table(note="STRINGdb/igraph not installed")
fwrite(hubs, opt$out_hubs, sep="\t")
message("[enrich] wrote ", opt$out_enrich, " and ", opt$out_hubs)
