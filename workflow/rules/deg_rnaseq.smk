# DESeq2 on pre-computed counts (RNA-seq strata).
rule deg_rnaseq:
    input: "data/expr/{dataset}.expr.rds"
    output: "results/deg/{dataset}.deg.tsv"
    params:
        group=lambda wc: samples.loc[wc.dataset, "group_column"],
        case=lambda wc: samples.loc[wc.dataset, "case_label"],
        ctrl=lambda wc: samples.loc[wc.dataset, "control_label"],
        padj=config["deg"]["padj"], lfc=config["deg"]["lfc"],
    conda: "../../envs/deseq2.yaml"
    resources: mem_mb=4000
    log: "results/logs/deg_rnaseq_{dataset}.log"
    shell:
        r"""
        Rscript workflow/scripts/deg_deseq2.R --expr {input} \
            --group {params.group} --case "{params.case}" --control "{params.ctrl}" \
            --padj {params.padj} --lfc {params.lfc} --out {output} > {log} 2>&1
        """
# NOTE: snakemake picks deg_rnaseq OR deg_microarray per dataset via the technique
#       column; wire the selection with a small input function if you prefer one rule.
