# Gut biopsies, PCoA coloured by site.
source("R/setup.R")
tses <- load_tses()

save_manuscript_figure(fig_beta_diversity_gut_site(build_beta_gut_site(tses$gut)),
                       "beta_diversity_gut_site")
