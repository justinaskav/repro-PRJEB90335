# Regenerate every manuscript figure into figures/. Each script is self-contained.
for (f in sort(list.files("R/figures", pattern = "\\.R$", full.names = TRUE))) {
  message("== ", f)
  source(f, local = new.env())
}
