# PRJEB90335 - reproduction code

Regenerates the manuscript figures and results from the raw reads at ENA
([PRJEB90335](https://www.ebi.ac.uk/ena/browser/view/PRJEB90335)).

```bash
git clone --recurse-submodules https://github.com/justinaskav/repro-PRJEB90335
cd repro-PRJEB90335

python3 scripts/download_data.py      # reads + GRCh38, checksum-verified
nextflow run workflows/main.nf        # merged table, taxonomy, tree
nextflow run workflows/r_prepare.nf   # TSE objects + diversity
nextflow run workflows/figures.nf     # figures/ (or: Rscript R/make_all.R)
```

Nextflow builds the conda environments from `envs/qiime2.lock.txt` and `envs/r.lock.txt` on the
first run; pass `--qiime2_conda_env <path> --r_conda_env <path>` to use existing ones.

To start from the processed data instead, put the TSE objects from Zenodo
([10.5281/zenodo.23318201](https://doi.org/10.5281/zenodo.23318201)) in `_output/r/tses/final/`
and run only the figure step. `Rscript scripts/export_processed_data.R` rebuilds that archive
from the outputs.

## Notes

- Three NovaSeq 6000 runs (`batch1`–`batch3`), denoised separately. 281 runs deposited, 279
  analysed (`data/metadata/excluded_samples.txt`); column descriptions in `data/metadata/README.md`.
- Sample IDs are `<participant>.<site>`: `IL` ileum, `PC` proximal colon, `DC` distal colon,
  `LE` lesion, `FE` feces. `Age` is the midpoint of a 5-year band.
- Reads carry the 8-nt barcode and primer as submitted; cutadapt trims them using `data/barcodes/`.
