suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(hdf5r)
  library(grid)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")

project_root <- normalizePath(getwd())
result_root <- file.path(project_root, "results")
h5ad_root <- file.path(
  project_root,
  "resources",
  "Brain_33558695_h5ad"
)

out_dir <- file.path(
  project_root,
  "plot",
  "human_dlpfc_cell_type_enrichment"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_order <- c("ADHD", "ANX", "BIP", "MDD", "SCZ")
samples <- paste0("V", 1:12)

sample_pairs <- tibble(
  sample = samples,
  pair = rep(paste0("Pair", 1:6), each = 2)
)

pair_order <- unique(sample_pairs$pair)

selected_gwas <- c(
  ADHD = "ADHD2022_iPSYCH_deCODE_PGC",
  ANX = "ANX_2026_daner_fullANX_v12_woUTAH_11022026",
  BIP = "bip2024_eur_noUKB_no23andMe",
  MDD = "pgc_mdd2025_no23andMe_noUKBB_eur_v3_49_24_11",
  SCZ = "PGC3_SCZ_wave3_european_autosome_public_v3_vcf"
)

sr_order <- c("SR7", "SR5", "SR6", "SR1", "SR4", "SR2", "SR3")

# Use the same symmetric log2(OR) color scale for both plots: [-3, 3] corresponds to OR [1/8, 8].
log2_or_limit <- 3

write_table <- function(x, name) {
  write_tsv(x, file.path(out_dir, name), na = "NA")
  invisible(x)
}

read_obs_column <- function(h5, column) {
  obs <- h5[["obs"]]
  x <- obs[[column]]

  if (inherits(x, "H5Group") && all(c("codes", "categories") %in% names(x))) {
    codes <- as.integer(x[["codes"]]$read())
    categories <- as.character(x[["categories"]]$read())
    values <- rep(NA_character_, length(codes))
    keep <- codes >= 0
    values[keep] <- categories[codes[keep] + 1L]
    return(values)
  }

  as.character(x$read())
}

read_spot_metadata <- function() {
  files <- file.path(
    h5ad_root,
    sprintf("HS_Brain_33558695_%s.seuratobj.h5ad", samples)
  )

  map_dfr(files, function(path) {
    h5 <- H5File$new(path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)

    sample_id <- str_match(
      basename(path),
      "_(V[0-9]+)\\.seuratobj\\.h5ad$"
    )[, 2]

    sr_cell_type <- read_obs_column(h5, "spatial_region_cell_type")

    tibble(
      sample = sample_id,
      barcode_raw = read_obs_column(h5, "id_raw"),
      sr_cell_type = sr_cell_type,
      sr = paste0("SR", str_extract(sr_cell_type, "^[^.]+")),
      cell_type = str_remove(sr_cell_type, "^[^.]+\\.")
    )
  }) %>%
    filter(
      !is.na(sr_cell_type),
      sr_cell_type != "nan",
      str_detect(sr_cell_type, "\\."),
      sr %in% sr_order,
      !is.na(cell_type),
      cell_type != "nan"
    ) %>%
    distinct(sample, barcode_raw, sr, cell_type, sr_cell_type) %>%
    arrange(sample, barcode_raw, match(sr, sr_order)) %>%
    group_by(sample, barcode_raw) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    left_join(sample_pairs, by = "sample")
}

parse_enrichment <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE, comment = "Note:")
  names(x) <- c(
    "Condition",
    "Median_IQR",
    "EnrichmentScore",
    "Adjusted_p"
  )

  x %>%
    mutate(
      EnrichmentScore = as.numeric(EnrichmentScore),
      Adjusted_p = as.numeric(Adjusted_p),
      barcode_raw = sub(
        "^([^.]+\\.[^.]+\\.[^.]+).*",
        "\\1",
        Condition
      )
    ) %>%
    group_by(barcode_raw) %>%
    slice_min(Adjusted_p, n = 1, with_ties = FALSE) %>%
    ungroup()
}

read_primary_spots <- function(spot_meta) {
  primary_catalog <- tibble(
    disease = disease_order,
    gwas = unname(selected_gwas[disease_order])
  )

  pmap_dfr(primary_catalog, function(disease, gwas) {
    map_dfr(samples, function(sample_id) {
      path <- file.path(
        result_root,
        disease,
        gwas,
        sample_id,
        sprintf(
          "%s_gene_marker_score.feather.enrichment.tsv",
          sample_id
        )
      )

      parse_enrichment(path) %>%
        mutate(
          disease = disease,
          gwas = gwas,
          sample = sample_id,
          .before = 1
        )
    })
  }) %>%
    left_join(sample_pairs, by = "sample") %>%
    mutate(is_enriched_spot = Adjusted_p < 0.05) %>%
    inner_join(
      spot_meta,
      by = c("sample", "barcode_raw", "pair")
    )
}

# Calculate the common Mantel-Haenszel OR stratified by pair.
descriptive_mh_or <- function(df, in_group) {
  strata <- intersect(pair_order, unique(as.character(df$pair)))

  components <- map_dfr(strata, function(pair_id) {
    keep <- as.character(df$pair) == pair_id
    enriched <- df$is_enriched_spot[keep]
    inside <- in_group[keep]

    a <- sum(enriched & inside)
    b <- sum(enriched & !inside)
    c <- sum(!enriched & inside)
    d <- sum(!enriched & !inside)
    n <- a + b + c + d

    tibble(
      numerator_component = a * d / n,
      denominator_component = b * c / n
    )
  })

  numerator <- sum(components$numerator_component)
  denominator <- sum(components$denominator_component)

  odds_ratio <- case_when(
    numerator == 0 && denominator == 0 ~ NA_real_,
    denominator == 0 ~ Inf,
    numerator == 0 ~ 0,
    TRUE ~ numerator / denominator
  )

  tibble(
    odds_ratio = odds_ratio,
    estimator = "Pair-stratified descriptive MH OR"
  )
}

summarise_enrichment <- function(df, group_column, groups) {
  map_dfr(disease_order, function(disease_id) {
    disease_data <- filter(df, disease == disease_id)

    map_dfr(groups, function(group_id) {
      inside <- as.character(disease_data[[group_column]]) == group_id
      enriched <- disease_data$is_enriched_spot

      descriptive_mh_or(disease_data, inside) %>%
        mutate(
          disease = disease_id,
          stratum = group_id,
          enriched_in = sum(enriched & inside),
          enriched_out = sum(enriched & !inside),
          nonenriched_in = sum(!enriched & inside),
          nonenriched_out = sum(!enriched & !inside),
          prop_enriched_in = enriched_in /
            max(enriched_in + nonenriched_in, 1),
          prop_enriched_out = enriched_out /
            max(enriched_out + nonenriched_out, 1)
        )
    })
  }) %>%
    mutate(
      log2_odds_ratio = log2(odds_ratio),
      log2_odds_ratio_plot = log2(
        pmin(
          pmax(odds_ratio, 2^(-log2_or_limit)),
          2^log2_or_limit
        )
      )
    )
}

# Perform enrichment analysis by cell_type.

spot_meta <- read_spot_metadata()
all_primary <- read_primary_spots(spot_meta)

# Order all cell types from highest to lowest by the number of spots.
cell_type_order <- all_primary %>%
  distinct(sample, barcode_raw, cell_type) %>%
  count(cell_type, sort = TRUE) %>%
  pull(cell_type)

cell_spots <- all_primary %>%
  distinct(
    disease,
    sample,
    pair,
    barcode_raw,
    cell_type,
    is_enriched_spot
  )

all_cell_data <- summarise_enrichment(
  cell_spots,
  group_column = "cell_type",
  groups = cell_type_order
) %>%
  mutate(
    disease = factor(disease, levels = disease_order),
    stratum = factor(stratum, levels = cell_type_order)
  )

write_table(
  all_cell_data,
  "human_dlpfc_cell_type_enrichment.tsv"
)

# Plot all cell types.

theme_heat <- theme_minimal(
  base_size = 22,
  base_family = "Arial"
) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(
      size = 25,
      color = "black",
      face = "bold"
    ),
    axis.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18, color = "black"),
    legend.key.height = unit(1.1, "cm"),
    legend.key.width = unit(0.65, "cm"),
    plot.margin = margin(3, 3, 3, 3)
  )

cell_group <- function(x) {
  case_when(
    str_detect(str_to_lower(x), "excit") ~ "Excitatory",
    str_detect(
      str_to_lower(x),
      "inhibit|pvalb|sst|sv2c|cxcl14"
    ) ~ "Inhibitory",
    TRUE ~ "Other"
  )
}

cell_label <- function(x) {
  case_when(
    x == "L2/3 excitatory neuron" ~ "L2/3 excitatory",
    x == "L4 excitatory neuron" ~ "L4 excitatory",
    x == "PLCH1 L4/5 excitatory neuron" ~ "PLCH1 L4/5 excitatory",
    x == "TSHZ2 L4/5 excitatory neuron" ~ "TSHZ2 L4/5 excitatory",
    x == "L5b excitatory neuron" ~ "L5b excitatory",
    x == "L5/6 excitatory neuron" ~ "L5/6 excitatory",
    x == "L6 excitatory neuron" ~ "L6 excitatory",
    x == "L6b excitatory neuron" ~ "L6b excitatory",
    x == "CXCL14 inhibitory neuron" ~ "CXCL14 inhibitor",
    x == "PVALB inhibitory neuron" ~ "PVALB inhibitory",
    x == "SST inhibitory neuron" ~ "SST inhibitory",
    x == "SV2C inhibitory neuron" ~ "SV2C inhibitory",
    x == "Oligodendrocyte progenitor cell" ~
      "Oligodendrocyte progenitor cell",
    x == "Oligodendrocyte" ~ "Oligodendrocyte",
    x == "Fibrous astrocyte" ~ "Fibrous astrocyte",
    x == "Pyramidal neuron" ~ "Pyramidal neuron",
    x == "Endothelial cell" ~ "Endothelial cell",
    x == "Microglia" ~ "Microglia",
    TRUE ~ str_replace_all(x, " neuron| cell", "")
  )
}

# Retain all cell types and sort them by the Excitatory, Inhibitory, and Other groups.
# Within each group, sort by the mean descriptive log2(OR) across the five diseases.
all_cell_types_ordered <- all_cell_data %>%
  mutate(cell_group = cell_group(as.character(stratum))) %>%
  group_by(stratum) %>%
  summarise(
    cell_group = first(cell_group),
    mean_log2_or = mean(log2_odds_ratio_plot, na.rm = TRUE),
    max_log2_or = max(log2_odds_ratio_plot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cell_group = factor(
      cell_group,
      levels = c("Excitatory", "Inhibitory", "Other")
    )
  ) %>%
  arrange(cell_group, desc(mean_log2_or), desc(max_log2_or)) %>%
  pull(stratum) %>%
  as.character()

plot_data <- all_cell_data %>%
  mutate(
    disease_plot = factor(
      as.character(disease),
      levels = rev(disease_order)
    ),
    stratum = factor(
      as.character(stratum),
      levels = all_cell_types_ordered
    )
  )

p <- ggplot(
  plot_data,
  aes(stratum, disease_plot, fill = log2_odds_ratio_plot)
) +
  geom_tile(color = "white", linewidth = 0.65) +
  scale_x_discrete(labels = cell_label, drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient2(
    low = "#263B73",
    mid = "#FFFFFF",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-log2_or_limit, log2_or_limit),
    breaks = seq(-log2_or_limit, log2_or_limit, by = 1),
    na.value = "#BDBDBD",
    name = "log2(OR)",
    guide = guide_colorbar(title.position = "top")
  ) +
  labs(x = "PSC cell type", y = NULL) +
  coord_fixed() +
  theme_heat +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1),
    axis.text.y = element_text(),
    axis.title.x = element_text()
  )

plot_width <- max(16, 4 + 0.65 * length(all_cell_types_ordered))

stem <- file.path(out_dir, "human_dlpfc_cell_type_enrichment")

ggsave(
  paste0(stem, ".png"),
  p,
  width = plot_width,
  height = 7.4,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  paste0(stem, ".pdf"),
  p,
  width = plot_width,
  height = 7.4,
  device = grDevices::cairo_pdf,
  bg = "white",
  limitsize = FALSE
)