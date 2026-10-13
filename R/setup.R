suppressPackageStartupMessages({
  library(tidyverse)
  library(mia)
  library(patchwork)
})

if (!exists("project_root")) project_root <- "."

shared <- file.path(project_root, "shared", "R")
source(file.path(shared, "util/classify_significance.R"))
source(file.path(shared, "util/save_figure.R"))
source(file.path(shared, "taxonomic/aggregate_taxa.R"))
source(file.path(shared, "taxonomic/resolve_ambiguous_taxa.R"))
source(file.path(shared, "taxonomic/order_taxa_by_abundance.R"))
source(file.path(shared, "taxonomic/calc_taxonomic_summary.R")) # calc_bootstrap_median_ci()
source(file.path(shared, "taxonomic/plot_composition_by_group.R"))
source(file.path(shared, "taxonomic/plot_dominance_bubble.R"))
source(file.path(shared, "alpha_diversity/calc_alpha_diversity_summary.R"))
source(file.path(shared, "alpha_diversity/alpha_plot_style.R"))
source(file.path(shared, "alpha_diversity/plot_alpha_diversity_summary_lmm.R"))
source(file.path(shared, "alpha_diversity/create_alpha_diversity_panel.R"))
source(file.path(shared, "alpha_diversity/plot_alpha_delta_vs_ref.R"))
source(file.path(shared, "alpha_diversity/plot_rarefaction_curves.R"))
source(file.path(shared, "beta_diversity/beta_diversity.R"))
source(file.path(shared, "beta_diversity/create_beta_diversity_panel.R"))
source(file.path(shared, "beta_diversity/combine_diversity_panels.R"))
source(file.path(shared, "differential_abundance/prepare_maaslin3_data.R"))
source(file.path(shared, "differential_abundance/plot_maaslin3_volcano.R"))
source(file.path(shared, "differential_abundance/drop_denominator_features.R"))
source(file.path(shared, "differential_abundance/plot_emm_interaction.R"))

source(file.path(project_root, "R", "helpers.R"))

theme_set(theme_project())
set.seed(2026) # same global seed as the report
