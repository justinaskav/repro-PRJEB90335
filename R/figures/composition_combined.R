# Mean composition per diagnosis x site: (a) phyla, (b) top 10 genera shaded by phylum.
source("R/setup.R")
tses <- load_tses()

phylum <- prep_phylum(tses)
genus <- prep_genus_all(tses)
save_manuscript_figure(
  fig_composition_combined(phylum$base_phylum, genus$base_genus_all, genus$tse_genus_all),
  "composition_combined"
)
