# Standardized DEG output schema

Every per-dataset DEG file MUST be a tab-separated file with exactly these columns,
so meta-analysis can read any of them the same way.

| column  | type   | meaning                                  |
|---------|--------|------------------------------------------|
| gene    | string | HGNC symbol (map probes/Ensembl first)   |
| log2FC  | float  | log2 fold change, case vs control        |
| p       | float  | raw p-value                              |
| padj    | float  | BH-adjusted p-value                      |
| n       | int    | total samples used in the contrast       |

Rules:
- Gene IDs harmonised to **HGNC symbols** before writing (one row per gene; collapse duplicate probes by max absolute log2FC).
- Direction is always **case minus control** (e.g., PCOS minus Control).
- File name: `results/deg/<dataset_id>.deg.tsv`.
- Do not pre-filter by significance here; meta step applies thresholds.
