# Cross-DATASET-validated gene panel (elastic net / RF) on meta-DEGs.
rule ml_panel:
    input: expand("results/meta/{stratum}.meta.tsv", stratum=STRATA_IDS)
    output: "results/ml/panel_report.tsv"
    conda: "../../envs/pyml.yaml"
    resources: mem_mb=4000
    log: "results/logs/ml_panel.log"
    shell:
        r"""
        python workflow/scripts/ml_panel.py --meta_dir results/meta \
            --expr_dir data/expr --samplesheet {config[samplesheet]} \
            --cv leave-one-dataset-out --out {output} > {log} 2>&1
        """
