suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(hdf5r)
  library(grid)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")
project_root <- normalizePath(getwd())
result_root <- file.path(project_root, "results")
h5ad_root <- file.path(project_root, "resources", "Brain_33558695_h5ad")
out_dir <- file.path(project_root, "plot", "human_dlpfc_profile_similarity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_order <- c("ADHD", "ANX", "BIP", "MDD", "SCZ")
samples <- paste0("V", 1:12)
sample_pairs <- tibble(
  sample = samples,
  pair = rep(paste0("Pair", 1:6), each = 2)
)

selected_gwas <- c(
  ADHD = "ADHD2022_iPSYCH_deCODE_PGC",
  ANX = "ANX_2026_daner_fullANX_v12_woUTAH_11022026",
  BIP = "bip2024_eur_noUKB_no23andMe",
  MDD = "pgc_mdd2025_no23andMe_noUKBB_eur_v3_49_24_11",
  SCZ = "PGC3_SCZ_wave3_european_autosome_public_v3_vcf"
)

sr_order <- c("SR7", "SR5", "SR6", "SR1", "SR4", "SR2", "SR3")

write_table <- function(x, name) {
  write_tsv(x, file.path(out_dir, name), na = "NA")
  invisible(x)
}

read_obs_column <- function(h5, column) {
  obs <- h5[["obs"]]
  x <- obs[[column]]

  if (
    inherits(x, "H5Group") &&
      all(c("codes", "categories") %in% names(x))
  ) {
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
    mutate(is_sig = Adjusted_p < 0.05) %>%
    inner_join(
      spot_meta,
      by = c("sample", "barcode_raw", "pair")
    )
}

spot_meta <- read_spot_metadata()
all_primary <- read_primary_spots(spot_meta)

spatial_profile <- all_primary %>%
  group_by(disease, sr) %>%
  summarise(
    mean_es = mean(EnrichmentScore, na.rm = TRUE),
    prop_sig = mean(is_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(feature = paste0("SR:", sr))

cell_profile <- all_primary %>%
  group_by(disease, cell_type) %>%
  summarise(
    mean_es = mean(EnrichmentScore, na.rm = TRUE),
    prop_sig = mean(is_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(feature = paste0("Cell:", cell_type))

profile_data <- bind_rows(
  spatial_profile %>% select(disease, feature, mean_es),
  cell_profile %>% select(disease, feature, mean_es)
)

profile_wide <- profile_data %>%
  pivot_wider(
    names_from = disease,
    values_from = mean_es,
    values_fill = 0
  )

correlation_matrix <- cor(
  as.matrix(profile_wide[, disease_order]),
  method = "spearman",
  use = "pairwise.complete.obs"
)

plot_data <- as.data.frame(as.table(correlation_matrix)) %>%
  rename(
    disease_1 = Var1,
    disease_2 = Var2,
    spearman_rho = Freq
  ) %>%
  mutate(
    disease_1 = factor(disease_1, levels = disease_order),
    disease_2 = factor(disease_2, levels = rev(disease_order)),
    label = sprintf("%.2f", spearman_rho)
  )

write_table(
  profile_data,
  "human_dlpfc_spatial_cellular_profiles.tsv"
)

write_table(
  plot_data,
  "human_dlpfc_profile_similarity.tsv"
)

theme_heat <- theme_minimal(
  base_size = 22,
  base_family = "Arial"
) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 22, color = "black"),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18, color = "black"),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5)
  )

plot_data <- plot_data %>%
  mutate(
    i = match(as.character(disease_1), disease_order),
    j = match(as.character(disease_2), disease_order)
  )

diagonal <- filter(plot_data, i == j)
lower_triangle <- filter(plot_data, i < j)

p <- ggplot() +
  geom_tile(
    data = lower_triangle,
    aes(disease_1, disease_2, fill = spearman_rho),
    color = "white",
    linewidth = 0.55
  ) +
  geom_tile(
    data = diagonal,
    aes(disease_1, disease_2),
    fill = "#F0F0F0",
    color = "white",
    linewidth = 0.55
  ) +
  geom_text(
    data = lower_triangle,
    aes(disease_1, disease_2, label = label),
    size = 7,
    color = "#171717"
  ) +
  geom_text(
    data = diagonal,
    aes(disease_1, disease_2, label = label),
    size = 7,
    color = "#171717"
  ) +
  scale_fill_gradient2(
    low = "#8FB6D8",
    mid = "white",
    high = "#B2182B",
    midpoint = 0.5,
    limits = c(0, 1),
    name = "Spearman rho",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barwidth = unit(0.8, "cm"),
      barheight = unit(3.8, "cm")
    )
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_heat +
  theme(
    legend.title = element_text(
      size = 18,
      face = "bold",
      margin = margin(b = 8)
    ),
    legend.text = element_text(
      size = 16,
      color = "black"
    ),

    legend.position = c(0.86, 0.73),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.margin = margin(3, 4, 3, 4),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 17
    ),
    axis.text.y = element_text(
      angle = 45,
      hjust = 1,
      size = 17
    )
  )

stem <- file.path(out_dir, "human_dlpfc_profile_similarity")

ggsave(
  paste0(stem, ".png"),
  p,
  width = 6.4,
  height = 5.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  paste0(stem, ".pdf"),
  p,
  width = 6.4,
  height = 5.8,
  device = grDevices::cairo_pdf,
  bg = "white"
)
