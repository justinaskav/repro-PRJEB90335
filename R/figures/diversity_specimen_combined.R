# Faeces vs mucosa: alpha (a Shannon, b Faith's PD) over beta (c Bray-Curtis, d unweighted UniFrac).
source("R/setup.R")
tses <- load_tses()

alpha_spec <- prep_alpha_specimen(tses)
beta_specimen <- build_beta_specimen(prep_tse_specimen(tses))
save_manuscript_figure(
  fig_diversity_specimen_combined(build_alpha_specimen_panel(alpha_spec$alpha_data_specimen),
                                  beta_specimen),
  "diversity_specimen_combined"
)
