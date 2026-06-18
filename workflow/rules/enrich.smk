# Functional enrichment + STRING PPI hub genes on the meta signature.
rule enrich:
    input: "results/meta/{stratum}.meta.tsv"
    output:
        enr="results/enrich/{stratum}.enrich.tsv",
        hubs="results/enrich/{stratum}.hubs.tsv"
    conda: "../../envs/meta.yaml"
    resources: mem_mb=3000
    log: "results/logs/enrich_{stratum}.log"
    shell:
        r"""
        Rscript workflow/scripts/enrich_hubs.R --meta {input} \
            --out_enrich {output.enr} --out_hubs {output.hubs} > {log} 2>&1
        """
