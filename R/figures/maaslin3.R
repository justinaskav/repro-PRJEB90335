# MaAsLin3 figures:
#   maaslin3_diagnosis_combined  (a) lesion vs ileum volcano, base model;
#                                (b) Fusobacterium EMMs, interaction model (all patients)
#   maaslin3_emm_supplementary   Parvimonas EMMs, same model
source("R/setup.R")
tses <- load_tses()

tse_full_genus <- maaslin3_genus_tse(tses$full)
maaslin3_full_genus <- run_maaslin3_genus(
  tse_full_genus, formula_str_full, "_output/r/maaslin3/full_genus"
)
maaslin3_genus_interaction <- run_maaslin3_genus(
  interaction_genus_tse(tse_full_genus), formula_str_interaction,
  "_output/r/maaslin3/genus_interaction"
)

save_manuscript_figure(
  fig_maaslin3_diagnosis_combined(
    build_maaslin_volcano_polyp(maaslin3_full_genus),
    build_emm_panel(maaslin3_genus_interaction, emm_main_taxon)
  ),
  "maaslin3_diagnosis_combined"
)

emm_supp_panels <- lapply(emm_supp_taxa, function(taxon) {
  build_emm_panel(maaslin3_genus_interaction, taxon)
})
save_manuscript_figure(fig_maaslin3_emm_supplementary(emm_supp_panels), "maaslin3_emm_supplementary")
