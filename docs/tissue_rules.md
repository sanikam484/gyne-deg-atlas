# Tissue / cell-type rules (read before adding any dataset)

These rules prevent the single most common failure in expression meta-analysis:
combining samples that are not biologically comparable.

## Hard rules
1. **One analysis per (condition x tissue) stratum.** Never merge tissues.
2. **Assign exactly ONE primary tissue per dataset** in the sample sheet (`tissue_primary`).
   The auto tissue tags overlap (a granulosa dataset may also say "ovarian"); pick the most specific.
3. **Human only.** Drop non-human / mixed-species series.
4. **RNA-seq and microarray analysed separately**, combined only at the results level.
5. **Case vs control must be the same tissue and same cell fraction.**

## Tissue distinctions that change the biology (do not blur)
- **Eutopic vs ectopic endometrium** (endometriosis/adenomyosis): eutopic = lining inside the uterus;
  ectopic = lesion. These are different comparisons. Record in `subtissue`. Never pool.
- **Granulosa vs cumulus vs theca vs whole ovary vs oocyte** (PCOS): related but distinct. Keep separate.
- **Endometrial epithelium vs stroma** (sorted) vs **whole endometrium** (bulk): keep separate.
- **Whole blood vs PBMC vs serum/plasma (cell-free):** different RNA pools. Keep separate.
- **Primary tissue vs cell line:** never mix.
- **Lesion vs adjacent normal vs healthy-control tissue:** define the contrast explicitly.

## Confounders to record in the sample sheet (so they can be modelled or matched)
- **Menstrual cycle phase** (proliferative vs secretory) - dominant for endometrium.
- **Ovarian stimulation / IVF** status - dominant for granulosa/cumulus.
- **BMI / obesity** - dominant for PCOS (esp. adipose, blood, muscle).
- **Age, medication (metformin, OCP), comorbidity.**

## Minimum-N rule
- A stratum is **meta-analysed only if it has >= `meta.min_datasets` (default 4)** datasets.
- Smaller strata -> single-dataset/descriptive DEG, reported as such (never called a "meta-analysis").
