# Gut biopsies by diagnosis: alpha (a Shannon, b Faith's PD) over beta (c Bray-Curtis, d unweighted UniFrac).
source("R/setup.R")
tses <- load_tses()

alpha <- prep_alpha(tses)
save_manuscript_figure(
  fig_diversity_diagnosis_combined(
    build_alpha_diagnosis_panel(alpha$alpha_data_gut, alpha$available_indices_gut),
    build_beta_gut_diagnosis(tses$gut)
  ),
  "diversity_diagnosis_combined"
)
