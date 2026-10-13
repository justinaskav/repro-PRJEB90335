
# --- Levels, palettes, shapes ------------------------------------------------------------

site_levels <- c("Ileum", "Proximal", "Distal", "Lesion", "Feces")

diagnosis_levels <- c("Control", "Adenoma", "Advanced adenoma", "Carcinoma")
diagnosis_colors <- c(
  "Control"          = "#1B9E77",
  "Adenoma"          = "#D95F02",
  "Advanced adenoma" = "#7570B3",
  "Carcinoma"        = "#E7298A"
)

site_colors <- c(
  "Ileum"    = "#4A8DC8",
  "Proximal" = "#3AA8A0",
  "Distal"   = "#3A8F5C",
  "Lesion"   = "#8B5CB8",
  "Feces"    = "#9A7553"
)
site_gut_colors <- site_colors[c("Ileum", "Proximal", "Distal", "Lesion")]

specimen_colors <- c(Mucosa = "#CC6677", Faeces = "#88CCEE")

taxonomy_palette <- c(
  "#77AADD",
  "#4EB265",
  "#EE8866",
  "#EEDD88",
  "#FFAABB",
  "#99DDFF",
  "#BBCC33",
  "#AAAA00",
  "#44BB99",
  "#FFCC66",
  "#DDCC77",
  "#CC6677",
  "#88CCEE",
  "#117733",
  "#AA4499",
  "#332288",
  "#999933",
  "#882255"
)

diagnosis_shapes <- c(
  "Control"          = 16,
  "Adenoma"          = 17,
  "Advanced adenoma" = 15,
  "Carcinoma"        = 18
)
site_shapes <- setNames(
  c(16, 17, 15, 18, 3, 8)[seq_along(names(site_gut_colors))],
  names(site_gut_colors)
)

# --- Model formulas ----------------------------------------------------------------------

formula_covariates_base <- "Diagnosis + Sex + Age_z"
formula_main_effects_gut <- "Diagnosis + Site + Sex + Age_z"
formula_str_gut <- paste("~", formula_main_effects_gut, "+ (1 | Patient)")
formula_fixef_feces <- paste("~", formula_covariates_base)
formula_str_full <- paste("~", formula_main_effects_gut, "+ (1 | Patient)")
formula_str_interaction <- paste("~ Diagnosis * Site + Sex + Age_z + (1 | Patient)")
formula_str_specimen <- "~ Specimen + Diagnosis + Sex + Age_z + (1 | Patient)"

# --- Figure style ------------------------------------------------------------------------

theme_project <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = ggplot2::rel(1.1)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1), color = "grey50"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.9, "lines"),
      legend.title = ggplot2::element_text(size = ggplot2::rel(1), face = "bold"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
      axis.title = ggplot2::element_text(size = ggplot2::rel(0.95), face = "bold"),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 7)),
      plot.margin = ggplot2::margin(10, 10, 10, 10, "pt"),
      panel.spacing = ggplot2::unit(0.6, "lines"),
      strip.text = ggplot2::element_text(size = ggplot2::rel(0.9), face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      plot.tag = ggplot2::element_text(size = ggplot2::rel(1.1), face = "bold"),
      plot.tag.position = c(0.01, 0.99)
    )
}

alpha_style <- alpha_plot_style()

theme_panel <- function(p) {
  p + cowplot::panel_border() + cowplot::background_grid() +
    ggplot2::theme(plot.margin = ggplot2::margin(10, 15, 10, 15))
}

FIGURES <- list(
  rarefaction_curves           = c(10, 5),
  composition_combined         = c(9.5, 8.5),
  dominance_bubble             = c(9, 6),
  diversity_specimen_combined  = c(9.5, 6),
  beta_diversity_gut_site      = c(10, 4.5),
  diversity_site_combined      = c(8, 3.6),
  diversity_diagnosis_combined = c(9.5, 6),
  maaslin3_diagnosis_combined  = c(9, 7),
  maaslin3_emm_supplementary   = c(6.5, 4)
)

save_manuscript_figure <- function(plot, name, dir = "figures") {
  size <- FIGURES[[name]]
  if (is.null(size)) stop("no size registered for figure '", name, "'")
  set.seed(2026)
  save_figure(relabel_for_print(plot), name, width = size[1], height = size[2], dir = dir)
}

print_labels <- c("Adenoma" = "Early adenoma", "Faeces" = "Feces")

print_factor <- function(x, levels = diagnosis_levels) {
  x <- factor(as.character(x), levels = levels)
  levels(x) <- recode_labels(levels(x), print_labels)
  x
}

relabel_for_print <- function(p, map = print_labels) {
  if (inherits(p, "patchwork")) {
    for (i in seq_along(p)) p[[i]] <- relabel_for_print(p[[i]], map)
    return(p)
  }
  if (!ggplot2::is_ggplot(p)) return(p)
  built <- ggplot2::ggplot_build(p)$plot$scales$scales
  p$scales <- p$scales$clone()
  for (b in built) {
    if (!inherits(b, "ScaleDiscrete")) next
    aes <- b$aesthetics[1]
    explicit <- which(vapply(p$scales$scales, function(s) aes %in% s$aesthetics, logical(1)))
    if (length(explicit) > 0) {
      s <- p$scales$scales[[explicit[1]]]$clone()
      s$labels <- relabeled_labels(s$labels, map)
      p$scales$scales[[explicit[1]]] <- s
    } else if (aes %in% c("x", "y", "colour", "fill", "shape", "linetype")) {
      default_scale <- utils::getFromNamespace(paste0("scale_", aes, "_discrete"), "ggplot2")
      p <- p + default_scale(labels = relabeled_labels(ggplot2::waiver(), map))
    }
  }
  p
}

recode_labels <- function(x, map) {
  if (!is.character(x)) return(x)
  hit <- x %in% names(map)
  x[hit] <- map[x[hit]]
  x
}

relabeled_labels <- function(old, map) {
  if (is.character(old)) return(recode_labels(old, map))
  if (is.function(old)) return(function(x) old(recode_labels(x, map)))
  if (inherits(old, "waiver")) return(function(x) recode_labels(x, map))
  old
}

# --- Data loading ------------------------------------------------------------------------

load_tses <- function(root = ".") {
  tse_dir <- file.path(root, "_output/r/tses/final")
  tses <- list(
    full  = readRDS(file.path(tse_dir, "tse_full_with_beta.rds")),
    gut   = readRDS(file.path(tse_dir, "tse_gut_with_beta.rds")),
    feces = readRDS(file.path(tse_dir, "tse_feces_with_beta.rds"))
  )
  for (nm in names(tses)) {
    SummarizedExperiment::colData(tses[[nm]])$Age_z <-
      as.numeric(scale(SummarizedExperiment::colData(tses[[nm]])$Age))
  }
  tses
}

specimen_of <- function(site) {
  factor(dplyr::if_else(site == "Feces", "Faeces", "Mucosa"), levels = c("Mucosa", "Faeces"))
}

# --- Taxonomic composition ---------------------------------------------------------------

prep_phylum <- function(tses) {
  tse_phylum_top8 <- aggregate_taxa(tses$full, rank = "phylum", top = 8, unclassified = "group")
  base_phylum <- mia::meltSE(tse_phylum_top8, assay.type = "counts",
                             add.col = c("Diagnosis", "Site", "Patient")) |>
    order_taxa_by_abundance(taxon_col = "FeatureID", abundance_col = "counts")
  list(tse_phylum_top8 = tse_phylum_top8, base_phylum = base_phylum)
}

prep_genus_all <- function(tses) {
  tse_genus_all <- aggregate_taxa(tses$full, rank = "genus", top = NULL,
                                  unclassified = "remove", incertae_sedis = "remove")
  base_genus_all <- mia::meltSE(tse_genus_all, assay.type = "counts",
                                add.col = c("Diagnosis", "Site", "Patient")) |>
    dplyr::group_by(SampleID) |>
    dplyr::mutate(rel_abundance = counts / sum(counts)) |>
    dplyr::ungroup()
  list(tse_genus_all = tse_genus_all, base_genus_all = base_genus_all)
}

comp_phyla <- c("Bacteroidota", "Bacillota", "Pseudomonadota", "Actinomycetota", "Fusobacteriota")
comp_phylum_cols <- c(
  Bacteroidota   = "#4A8DC8",
  Bacillota      = "#E08A2E",
  Pseudomonadota = "#1F8A70",
  Actinomycetota = "#7B5EA7",
  Fusobacteriota = "#C0143C",
  "Other phyla"  = "#D9D9D9"
)
comp_shades <- list(
  Bacteroidota   = c("#1D4E89", "#4A8DC8", "#A9CCEB"),
  Bacillota      = c("#A0520C", "#E08A2E", "#F5C98B"),
  Pseudomonadota = c("#0B4A3B", "#1F8A70", "#9FD6C3"),
  Fusobacteriota = c("#C0143C")
)
comp_n_genera <- 10
site_abbrev <- c(Ileum = "Ileum", Proximal = "Prox", Distal = "Dist", Lesion = "Lesion", Feces = "Feces")

fig_composition_combined <- function(base_phylum, base_genus_all, tse_genus_all) {
  sample_totals <- base_phylum |>
    dplyr::group_by(SampleID) |>
    dplyr::summarise(total = sum(counts), .groups = "drop")
  sample_meta <- dplyr::distinct(base_phylum, SampleID, Diagnosis, Site)
  comp_mean <- function(df, taxon_col) {
    calc_group_mean_composition(df, taxon_col, totals = sample_totals,
                                facet_var = "Diagnosis", x_var = "Site",
                                facet_levels = diagnosis_levels, x_levels = site_levels,
                                x_abbrev = site_abbrev, sample_meta = sample_meta) |>
      dplyr::mutate(x = factor(site_abbrev[as.character(Site)], levels = site_abbrev),
                    Diagnosis = print_factor(Diagnosis))
  }

  phylum_levels <- c(comp_phyla, "Other phyla")
  comp_phylum <- base_phylum |>
    dplyr::mutate(phylum = ifelse(as.character(FeatureID) %in% comp_phyla,
                                  as.character(FeatureID), "Other phyla")) |>
    comp_mean("phylum") |>
    dplyr::mutate(taxon = factor(taxon, levels = rev(phylum_levels)))

  genus_phylum <- as.data.frame(SummarizedExperiment::rowData(tse_genus_all)) |>
    tibble::rownames_to_column("FeatureID") |>
    dplyr::select(FeatureID, phylum)
  top_genera <- base_genus_all |>
    dplyr::left_join(sample_totals, by = "SampleID") |>
    dplyr::group_by(FeatureID) |>
    dplyr::summarise(m = mean(counts / total), .groups = "drop") |>
    dplyr::slice_max(m, n = comp_n_genera) |>
    dplyr::mutate(FeatureID = as.character(FeatureID)) |>
    dplyr::left_join(genus_phylum, by = "FeatureID") |>
    dplyr::mutate(phylum = factor(phylum, levels = comp_phyla)) |>
    dplyr::arrange(phylum, dplyr::desc(m))

  genus_cols <- shade_by_parent(top_genera$FeatureID, top_genera$phylum, comp_shades)
  genus_cols["Other genera"] <- comp_phylum_cols[["Other phyla"]]
  genus_levels <- names(genus_cols)

  comp_genus <- base_genus_all |>
    dplyr::mutate(genus = ifelse(as.character(FeatureID) %in% top_genera$FeatureID,
                                 as.character(FeatureID), "Other genera")) |>
    dplyr::select(SampleID, Diagnosis, Site, genus, counts) |>
    dplyr::bind_rows(
      base_genus_all |>
        dplyr::group_by(SampleID, Diagnosis, Site) |>
        dplyr::summarise(counts = sum(counts), .groups = "drop") |>
        dplyr::left_join(sample_totals, by = "SampleID") |>
        dplyr::transmute(SampleID, Diagnosis, Site, genus = "Other genera", counts = total - counts)
    ) |>
    comp_mean("genus") |>
    dplyr::mutate(taxon = factor(taxon, levels = rev(genus_levels)))

  comp_style <- theme_project() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title.position = "left",
      legend.text = ggtext::element_markdown(size = ggplot2::rel(0.9))
    )
  comp_n <- list(
    ggplot2::geom_text(data = dplyr::distinct(comp_phylum, Diagnosis, x, n),
                       ggplot2::aes(x = x, y = 1.015, label = n), inherit.aes = FALSE,
                       vjust = 0, size = 2.6, colour = "grey40"),
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                breaks = seq(0, 1, 0.25),
                                expand = ggplot2::expansion(mult = c(0, 0.07)))
  )
  p_a <- plot_group_composition(comp_phylum, comp_phylum_cols, phylum_levels, "Phylum",
                                legend_rows = 2, facet_var = "Diagnosis") + comp_style + comp_n
  p_b <- plot_group_composition(comp_genus, genus_cols, genus_levels, "Genus",
                                legend_rows = 3, facet_var = "Diagnosis") + comp_style + comp_n
  p_a / p_b + patchwork::plot_annotation(tag_levels = "a")
}

# --- Dominance ---------------------------------------------------------------------------

prep_dominance <- function(base_genus_all, threshold = 0.5) {
  berger_parker <- calc_berger_parker(base_genus_all,
                                      group_cols = c("SampleID", "Patient", "Site", "Diagnosis"))
  dominated_patients <- berger_parker |>
    dplyr::filter(berger_parker > threshold) |>
    dplyr::distinct(Patient) |>
    dplyr::pull(Patient)
  bp_bubble <- berger_parker |> dplyr::filter(Patient %in% dominated_patients)
  dom_genus_levels <- bp_bubble |>
    dplyr::count(dominant_genus, sort = TRUE) |>
    dplyr::pull(dominant_genus)
  list(
    berger_parker = berger_parker,
    dom_threshold = threshold,
    dominated_patients = dominated_patients,
    bp_bubble = bp_bubble,
    dom_genus_levels = dom_genus_levels,
    dominance_palette = palette_for_levels(dom_genus_levels, taxonomy_palette),
    n_dominated_samples = sum(berger_parker$berger_parker > threshold),
    n_dominated_patients = length(dominated_patients)
  )
}

fig_dominance_bubble <- function(dom) {
  strip <- function(x) gsub(" ", "\n", as.character(print_factor(x)))
  bp <- dplyr::mutate(dom$bp_bubble, Diagnosis = strip(Diagnosis))
  plot_dominance_bubble(bp, dom$dominance_palette, dom$dom_genus_levels,
                        site_levels = site_levels,
                        facet_levels = strip(diagnosis_levels),
                        threshold = dom$dom_threshold,
                        size_range = c(2, 5.5), subthreshold_size = 1.6) +
    ggplot2::labs(caption = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 3.5, shape = 21)))
}

# --- Alpha diversity ---------------------------------------------------------------------

alpha_indices <- c("Shannon", "Faith_PD")

prep_alpha <- function(tses) {
  alpha_data_full <- as.data.frame(SummarizedExperiment::colData(tses$full))
  alpha_data_gut <- as.data.frame(SummarizedExperiment::colData(tses$gut))
  alpha_data_feces <- as.data.frame(SummarizedExperiment::colData(tses$feces))

  alpha_configs <- list(
    full = list(data = alpha_data_full, formula = formula_str_full,
                groups = list(Diagnosis = list(ref = "Control"), Site = list(ref = "Ileum"))),
    gut = list(data = alpha_data_gut, formula = formula_str_gut,
               groups = list(Diagnosis = list(ref = "Control"), Site = list(ref = "Ileum"))),
    feces = list(data = alpha_data_feces, formula = formula_fixef_feces,
                 groups = list(Diagnosis = list(ref = "Control")))
  )

  alpha_results <- list()
  for (ds in names(alpha_configs)) {
    cfg <- alpha_configs[[ds]]
    alpha_results[[ds]] <- list()
    for (idx in alpha_indices) {
      alpha_results[[ds]][[idx]] <- list()
      for (grp in names(cfg$groups)) {
        alpha_results[[ds]][[idx]][[grp]] <- calc_alpha_diversity_summary(
          data = cfg$data, y_var = idx, group_var = grp,
          lmm_formula_str = paste0(idx, " ", cfg$formula),
          ref_group = cfg$groups[[grp]]$ref,
          post_hoc_method = "BH", digits = 2
        )
      }
    }
  }

  list(
    alpha_data_full = alpha_data_full, alpha_data_gut = alpha_data_gut,
    alpha_data_feces = alpha_data_feces,
    available_indices_full = intersect(alpha_indices, colnames(alpha_data_full)),
    available_indices_gut = intersect(alpha_indices, colnames(alpha_data_gut)),
    available_indices_feces = intersect(alpha_indices, colnames(alpha_data_feces)),
    alpha_configs = alpha_configs, alpha_results = alpha_results
  )
}

prep_alpha_specimen <- function(tses) {
  alpha_data_specimen <- as.data.frame(SummarizedExperiment::colData(tses$full)) |>
    dplyr::mutate(Specimen = specimen_of(Site))
  alpha_results_specimen <- list()
  for (idx in alpha_indices) {
    alpha_results_specimen[[idx]] <- calc_alpha_diversity_summary(
      data = alpha_data_specimen, y_var = idx, group_var = "Specimen",
      lmm_formula_str = paste0(idx, " ", formula_str_specimen),
      ref_group = "Mucosa", post_hoc_method = "BH", digits = 2
    )
  }
  list(alpha_data_specimen = alpha_data_specimen, alpha_results_specimen = alpha_results_specimen)
}

build_alpha_diagnosis_panel <- function(alpha_data_gut, available_indices_gut) {
  create_alpha_diversity_panel(
    data = alpha_data_gut,
    x_var = "Diagnosis",
    available_indices = available_indices_gut,
    lmm_formula_str = paste0("INDEX ", formula_str_gut),
    ref_group = "Control",
    x_label = "Diagnosis",
    color_values = diagnosis_colors,
    shape_var = "Site",
    show_lmm_label = FALSE,
    style = alpha_style,
    post_process = function(p) {
      theme_panel(p) + ggplot2::scale_x_discrete(labels = scales::label_wrap(10))
    }
  )
}

build_alpha_specimen_panel <- function(alpha_data_specimen) {
  create_alpha_diversity_panel(
    data = alpha_data_specimen,
    x_var = "Specimen",
    available_indices = intersect(alpha_indices, colnames(alpha_data_specimen)),
    lmm_formula_str = paste0("INDEX ", formula_str_specimen),
    ref_group = "Mucosa",
    x_label = "Specimen",
    color_values = specimen_colors,
    shape_var = "Diagnosis",
    show_lmm_label = FALSE,
    style = alpha_style,
    post_process = function(p) {
      if (identical(p$labels$y, "Faith_PD")) p <- p + ggplot2::labs(y = "Faith’s PD")
      theme_panel(p)
    }
  )
}

build_alpha_site_delta <- function(alpha_data_gut, alpha_results, available_indices_gut) {
  panels <- lapply(available_indices_gut, function(idx) {
    plot_alpha_delta_vs_ref(
      data = alpha_data_gut,
      index = idx,
      posthoc = alpha_results$gut[[idx]]$Site$posthoc,
      ref_site = "Ileum",
      sites = c("Proximal", "Distal", "Lesion"),
      colors = site_gut_colors,
      style = alpha_style
    ) |> theme_panel()
  })
  patchwork::wrap_plots(panels, ncol = 2) + patchwork::plot_annotation(tag_levels = "a")
}

fig_diversity_site_combined <- function(alpha_site_delta_tagged) alpha_site_delta_tagged

# --- Beta diversity ----------------------------------------------------------------------

prep_tse_specimen <- function(tses) {
  tse <- tses$full
  SummarizedExperiment::colData(tse)$Specimen <- specimen_of(SummarizedExperiment::colData(tse)$Site)
  tse
}

build_beta_panel <- function(tse, group_var, formula, color_values, shape_by) {
  create_beta_diversity_panel(
    tse = tse,
    group_var = group_var,
    formula = formula,
    strata_var = "Patient",
    color_by = group_var,
    shape_by = shape_by,
    permutations = 9999,
    metric_filter = c("bray", "unifrac_unweighted"),
    post_process = function(p) {
      theme_panel(p) +
        ggplot2::scale_color_manual(values = color_values) +
        ggplot2::scale_fill_manual(values = color_values)
    }
  )
}

build_beta_specimen <- function(tse_specimen) {
  build_beta_panel(tse_specimen, "Specimen", "Specimen + Diagnosis + Sex + Age_z",
                   specimen_colors, shape_by = "Diagnosis")
}
build_beta_gut_diagnosis <- function(tse_gut) {
  build_beta_panel(tse_gut, "Diagnosis", formula_main_effects_gut, diagnosis_colors, shape_by = "Site")
}
build_beta_gut_site <- function(tse_gut) {
  build_beta_panel(tse_gut, "Site", formula_main_effects_gut, site_gut_colors, shape_by = "Diagnosis")
}

fig_beta_diversity_gut_site <- function(beta_gut_site) wrap_beta_panel(beta_gut_site)

fig_diversity_specimen_combined <- function(alpha_specimen_panel, beta_specimen) {
  combine_alpha_beta_panels(alpha_specimen_panel, beta_specimen,
                            color_values = specimen_colors, shape_values = diagnosis_shapes,
                            color_title = "Specimen", shape_title = "Diagnosis")
}

fig_diversity_diagnosis_combined <- function(alpha_diagnosis_panel, beta_gut_diagnosis) {
  combine_alpha_beta_panels(alpha_diagnosis_panel, beta_gut_diagnosis,
                            color_values = diagnosis_colors, shape_values = site_shapes,
                            color_title = "Diagnosis", shape_title = "Site")
}

# --- Differential abundance (MaAsLin3) ---------------------------------------------------

maaslin3_genus_tse <- function(tse) {
  resolve_ambiguous_taxa(tse, rank = "genus") |>
    aggregate_taxa(rank = "genus", top = NULL, unclassified = "group", incertae_sedis = "remainder")
}

run_maaslin3_genus <- function(tse_genus, formula, output, seed = 2026) {
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  set.seed(seed)
  res <- maaslin3::maaslin3(
    input_data = tse_genus,
    output = output,
    formula = formula,
    min_prevalence = 0.1,
    max_significance = 0.05,
    small_random_effects = TRUE
  )
  drop_denominator_features(res)
}

interaction_genus_tse <- function(tse_full_genus) {
  tse <- tse_full_genus
  tse$Diagnosis <- relevel(factor(tse$Diagnosis), ref = "Adenoma")
  tse
}

build_maaslin_volcano_polyp <- function(maaslin3_full_genus) {
  n_sig_polyp <- prepare_maaslin3_data(
    maaslin3_full_genus, coefficient_names = "SiteLesion", model_type = "joint"
  ) |> dplyr::filter(qval < 0.05) |> nrow()
  plot_maaslin3_volcano(
    maaslin3_full_genus,
    coefficient_name = "SiteLesion",
    model_type = "joint",
    q_threshold = 0.05,
    max_labels = n_sig_polyp,
    title = "Lesion vs ileum"
  ) +
    ggplot2::scale_color_manual(
      values = c("Significant Up" = "#e41a1c", "Significant Down" = "#377eb8",
                 "Not Significant" = "grey50"),
      labels = c("Significant Up" = "Enriched at lesion (q < 0.05)",
                 "Significant Down" = "Depleted at lesion (q < 0.05)",
                 "Not Significant" = "Not significant"),
      name = NULL
    ) +
    ggplot2::labs(subtitle = NULL, x = "β (lesion vs ileum)",
                  y = expression(-log[10] ~ italic(q)))
}

emm_main_taxon <- "Fusobacterium"
emm_supp_taxa  <- c("Parvimonas")

build_emm_panel <- function(maaslin3_interaction, taxon) {
  fit <- maaslin3_interaction$fit_data_abundance$fits[[taxon]]
  if (is.null(fit)) return(NULL)
  safe_render(plot_emm_interaction(
    fit, taxon,
    group_var = "Diagnosis", x_var = "Site",
    group_levels = diagnosis_levels, x_levels = site_levels,
    colors = diagnosis_colors, shapes = diagnosis_shapes,
    separate_before = "Feces", x_label = "Site",
    y_label = expression(Estimated ~ log[2] ~ relative ~ abundance)
  ))
}

fig_maaslin3_diagnosis_combined <- function(maaslin_volcano_polyp, emm_main_panel) {
  patchwork::wrap_plots(list(maaslin_volcano_polyp, emm_main_panel), ncol = 1, heights = c(1, 1)) +
    patchwork::plot_annotation(tag_levels = "a")
}

fig_maaslin3_emm_supplementary <- function(emm_supp_panels) {
  tags <- if (length(emm_supp_panels) > 1) "a" else NULL
  patchwork::wrap_plots(emm_supp_panels, ncol = min(2, length(emm_supp_panels)), guides = "collect") +
    patchwork::plot_annotation(tag_levels = tags) &
    ggplot2::theme(legend.position = "right")
}

# --- Rarefaction -------------------------------------------------------------------------

fig_rarefaction_curves <- function(tse_full) {
  calc_rarefaction_curves(tse_full, meta_cols = c("Diagnosis", "Site")) |>
    plot_rarefaction_curves(color_by = "Diagnosis", colors = diagnosis_colors)
}
