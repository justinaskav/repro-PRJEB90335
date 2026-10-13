# Build the Zenodo archive of processed data from this repo's outputs:
#   _output/processed_data/ and _output/processed_data.zip
# Run after r_prepare.nf and the MaAsLin3 figure (R/figures/maaslin3.R):
#   Rscript scripts/export_processed_data.R
suppressPackageStartupMessages({
  library(mia)
  library(Biostrings)
  library(ape)
})

tse_dir <- "_output/r/tses/final"
maaslin_dir <- "_output/r/maaslin3"
out <- "_output/processed_data"

# Published sample columns; everything else in colData (sequencing depths, derived groupings)
# stays out of the archive.
sample_cols <- c("SampleID", "Group", "Patient", "Age", "Age_band", "Sex", "DNA_concentration",
                 "Merge_Batch", "Site", "Location", "Diagnosis")
diversity_cols <- c("Shannon", "Faith_PD")

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

unlink(out, recursive = TRUE)
dir.create(file.path(out, "results"), recursive = TRUE)

for (subset in c("full", "gut", "feces")) {
  name <- paste0("tse_", subset, "_with_beta.rds")
  tse <- readRDS(file.path(tse_dir, name))
  missing <- setdiff(c(sample_cols, diversity_cols), names(colData(tse)))
  if (length(missing)) stop(name, " lacks colData columns: ", paste(missing, collapse = ", "))
  colData(tse) <- colData(tse)[, c(sample_cols, diversity_cols)]
  saveRDS(tse, file.path(out, name))
  if (subset == "full") full <- tse
}

asvs <- sort(rownames(full))
samples <- sort(colnames(full))

counts <- assay(full, "counts")[asvs, samples]
write_tsv(data.frame(ASV = asvs, counts, check.names = FALSE), file.path(out, "asv_counts.tsv"))

taxonomy <- as.data.frame(rowData(full)[asvs, ])
write_tsv(data.frame(ASV = asvs, taxonomy, check.names = FALSE), file.path(out, "asv_taxonomy.tsv"))

writeXStringSet(referenceSeq(full)[asvs], file.path(out, "asv_sequences.fasta"), width = 80)
write.tree(rowTree(full), file.path(out, "asv_tree.nwk"))

meta <- as.data.frame(colData(full)[samples, setdiff(sample_cols, "Group")])
write_tsv(meta, file.path(out, "sample_metadata.tsv"))

# MaAsLin3 output as written, including the denominator-only placeholder features
results <- c(full_genus = "maaslin3_main_effects_all_results.tsv",
             genus_interaction = "maaslin3_diagnosis_x_site_all_results.tsv")
for (run in names(results)) {
  src <- file.path(maaslin_dir, run, "all_results.tsv")
  if (!file.exists(src)) stop(src, " not found; run R/figures/maaslin3.R first")
  file.copy(src, file.path(out, "results", results[[run]]), overwrite = TRUE)
}

invisible(file.copy("scripts/processed_data_README.md", file.path(out, "README.md"), overwrite = TRUE))

zip_path <- file.path(normalizePath(dirname(out)), paste0(basename(out), ".zip"))
unlink(zip_path)
withr::with_dir(dirname(out), utils::zip(zip_path, basename(out), flags = "-rqX", zip = "zip"))
message("wrote ", out, "/ and ", zip_path)
