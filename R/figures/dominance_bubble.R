# Patient x site dominance (Berger-Parker), patients with at least one dominated sample.
source("R/setup.R")
tses <- load_tses()

genus <- prep_genus_all(tses)
save_manuscript_figure(fig_dominance_bubble(prep_dominance(genus$base_genus_all)), "dominance_bubble")
