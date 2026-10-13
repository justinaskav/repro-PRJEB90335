# Rarefaction curves, coloured by diagnosis.
source("R/setup.R")
tses <- load_tses()

save_manuscript_figure(fig_rarefaction_curves(tses$full), "rarefaction_curves")
