# Sample metadata

| File | Used by | Contents |
|---|---|---|
| `metadata.tsv` | `r_prepare.nf` | one row per sample, keyed by `SampleID`: `Group`, `Patient`, `Gut_Site`, `Age`, `Age_band`, `Sex`, `DNA_concentration`, `Merge_Batch` |
| `batch1.tsv`, `batch2.tsv`, `batch3.tsv` | `main.nf` | `SampleID` and `DNA_concentration` per sequencing run: which samples each run keeps, and the input of the decontam frequency method |
| `excluded_samples.txt` | `r_prepare.nf` | deposited samples left out of the analysis |

- `Group`: `K` control, `A1` early adenoma, `A2` advanced adenoma, `C` carcinoma (in situ or invasive).
- `Gut_Site`: `Ileum`, `Proximal`, `Distal`, `Lesion`, `Excrement` (feces).
- `Age`: midpoint of the 5-year band in `Age_band`. Exact ages are not published.
- `Merge_Batch`: sequencing run (1–3 = `batch1`–`batch3`).

Any further column in `metadata.tsv` is carried into the TSE `colData` unchanged.
