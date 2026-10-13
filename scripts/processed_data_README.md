# Processed data: lesion-confined *Fusobacterium* enrichment in the colorectal mucosa

Processed 16S rRNA (V3–V4) data for the study "Lesion-confined Fusobacterium enrichment and no
detectable field effect in the colorectal mucosa: a within-patient multi-site biopsy study".
Raw reads: ENA PRJEB90335. Code: https://github.com/justinaskav/repro-PRJEB90335, which also
builds this archive (`scripts/export_processed_data.R`).

279 samples from 60 participants (mucosal biopsies from the terminal ileum, proximal colon,
distal colon and lesion, plus stool); 12,216 amplicon sequence variants (ASVs).

## Files

| File | Content |
|---|---|
| `asv_counts.tsv` | ASV × sample read counts (DADA2, after host and contaminant removal) |
| `asv_taxonomy.tsv` | Taxonomy per ASV (SILVA 138.2, Naive Bayes classifier, V3–V4) |
| `asv_sequences.fasta` | ASV sequences |
| `asv_tree.nwk` | Rooted phylogenetic tree of the ASVs (FastTree2) |
| `sample_metadata.tsv` | De-identified sample metadata (see below) |
| `tse_full_with_beta.rds`, `tse_gut_with_beta.rds`, `tse_feces_with_beta.rds` | The same data as TreeSummarizedExperiment objects (R/Bioconductor): all samples, mucosal biopsies, stool, with alpha diversity (`Shannon`, `Faith_PD`) and beta-diversity ordinations. These are the objects the published results were computed from; place them in `_output/r/tses/final/` of the code repository to regenerate the figures. |
| `results/maaslin3_main_effects_all_results.tsv` | MaAsLin 3 results for all genera, model `~ Diagnosis + Site + Sex + Age_z + (1 | Patient)` |
| `results/maaslin3_diagnosis_x_site_all_results.tsv` | MaAsLin 3 results for all genera, model `~ Diagnosis * Site + Sex + Age_z + (1 | Patient)`, reference early adenoma × ileum |

## Sample metadata

- `Patient`: pseudonymous participant code; `Site`: Ileum, Proximal, Distal, Lesion or Feces.
- `Diagnosis`: Control, Adenoma (early adenoma), Advanced adenoma, Carcinoma; `Group` is the same as
  a code (K, A1, A2, C).
- `Age`: midpoint of a 5-year age band (`Age_band`); exact ages are not published. `Sex`: M/F.
- `DNA_concentration`: post-extraction DNA concentration (ng/µl), used for decontam.
- `Merge_Batch`: sequencing run (1–3). `Location`: Gut (mucosal biopsy) or Feces.

No questionnaire, clinical or immunohistochemical data are included.

## Licence

CC BY 4.0.
