# Per-patient change in alpha diversity relative to the patient's own ileum.
source("R/setup.R")
tses <- load_tses()

alpha <- prep_alpha(tses)
save_manuscript_figure(
  fig_diversity_site_combined(
    build_alpha_site_delta(alpha$alpha_data_gut, alpha$alpha_results, alpha$available_indices_gut)
  ),
  "diversity_site_combined"
)
