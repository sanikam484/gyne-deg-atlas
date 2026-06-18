# Contributing (subordinate workflow)

1. Pick your assigned (condition x tissue) stratum from `config/strata.csv`.
2. Fill your datasets into `config/samplesheet.csv` (copy from the template). Set:
   condition, tissue_primary, subtissue, technique, data_source, group_column, case_label,
   control_label, cycle_phase, notes. Confirm each against the GEO record.
3. Run a SINGLE-dataset smoke test first:
   `snakemake --use-conda --cores 1 results/deg/<your_first_dataset>.deg.tsv`
4. Commit your DEG tables + a short note in the PR describing QC (PCA, outliers, group sizes).
5. One branch per person: `feature/<condition>-<tissue>-<yourname>`. Open a PR; do not push to main.
6. Every commit is your visible GitHub contribution - commit small and often.
