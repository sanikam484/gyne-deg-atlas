# limma on GEO processed expression (microarray strata).
rule deg_microarray:
    input: "data/expr/{dataset}.expr.rds"
    output: "results/deg/{dataset}.deg.microarray.tsv"
    params:
        group=lambda wc: samples.loc[wc.dataset, "group_column"],
        case=lambda wc: samples.loc[wc.dataset, "case_label"],
        ctrl=lambda wc: samples.loc[wc.dataset, "control_label"],
    conda: "../../envs/limma.yaml"
    resources: mem_mb=3000
    log: "results/logs/deg_micro_{dataset}.log"
    shell:
        r"""
        Rscript workflow/scripts/deg_limma.R --expr {input} \
            --group {params.group} --case "{params.case}" --control "{params.ctrl}" \
            --out {output} > {log} 2>&1
        """
