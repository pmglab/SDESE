suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(hdf5r)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")
project_root <- normalizePath(getwd())
result_root <- file.path(project_root, "results")
h5ad_root <- file.path(project_root, "resources", "Brain_33558695_h5ad")
out_dir <- file.path(project_root, "plot", "human_dlpfc_gwas_profile_consistency")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_order <- c("ADHD", "ANX", "BIP", "MDD", "SCZ")
samples <- paste0("V", 1:12)
sample_pairs <- tibble(sample = samples, pair = rep(paste0("Pair", 1:6), each = 2))
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
  files <- file.path(h5ad_root, sprintf("HS_Brain_33558695_%s.seuratobj.h5ad", samples))
  map_dfr(files, function(path) {
    h5 <- H5File$new(path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    sample_id <- str_match(basename(path), "_(V[0-9]+)\\.seuratobj\\.h5ad$")[, 2]
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
      !is.na(sr_cell_type), sr_cell_type != "nan", str_detect(sr_cell_type, "\\."),
      sr %in% sr_order, !is.na(cell_type), cell_type != "nan"
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
  names(x) <- c("Condition", "Median_IQR", "EnrichmentScore", "Adjusted_p")
  x %>%
    mutate(
      EnrichmentScore = as.numeric(EnrichmentScore),
      Adjusted_p = as.numeric(Adjusted_p),
      barcode_raw = sub("^([^.]+\\.[^.]+\\.[^.]+).*", "\\1", Condition)
    ) %>%
    group_by(barcode_raw) %>%
    slice_min(Adjusted_p, n = 1, with_ties = FALSE) %>%
    ungroup()
}

spot_meta <- read_spot_metadata()
gwas_catalog <- map_dfr(disease_order, function(disease_id) {
  disease_dir <- file.path(result_root, disease_id)
  tibble(disease = disease_id, gwas = list.dirs(disease_dir, recursive = FALSE, full.names = FALSE))
}) %>%
  filter(gwas != "") %>%
  mutate(is_primary = gwas == selected_gwas[disease])

all_gwas <- pmap_dfr(gwas_catalog %>% select(disease, gwas), function(disease, gwas) {
  map_dfr(samples, function(sample_id) {
    path <- file.path(
      result_root, disease, gwas, sample_id,
      sprintf("%s_gene_marker_score.feather.enrichment.tsv", sample_id)
    )
    parse_enrichment(path) %>%
      mutate(disease = disease, gwas = gwas, sample = sample_id, .before = 1)
  })
}) %>%
  left_join(sample_pairs, by = "sample") %>%
  mutate(is_sig = Adjusted_p < 0.05) %>%
  inner_join(spot_meta, by = c("sample", "barcode_raw", "pair"))

spatial_profiles <- all_gwas %>%
  group_by(disease, gwas, sr) %>%
  summarise(
    mean_es = mean(EnrichmentScore, na.rm = TRUE),
    prop_sig = mean(is_sig, na.rm = TRUE),
    .groups = "drop"
  )
cell_profiles <- all_gwas %>%
  group_by(disease, gwas, cell_type) %>%
  summarise(
    mean_es = mean(EnrichmentScore, na.rm = TRUE),
    prop_sig = mean(is_sig, na.rm = TRUE),
    .groups = "drop"
  )

profile_consistency <- function(profile_data, feature_column) {
  map_dfr(disease_order, function(disease_id) {
    primary_gwas <- selected_gwas[[disease_id]]
    reference <- profile_data %>%
      filter(disease == disease_id, gwas == primary_gwas) %>%
      select(feature = all_of(feature_column), ref_es = mean_es)
    profile_data %>%
      filter(disease == disease_id) %>%
      select(gwas, feature = all_of(feature_column), mean_es) %>%
      inner_join(reference, by = "feature") %>%
      group_by(gwas) %>%
      summarise(
        spearman_rho = suppressWarnings(cor(
          mean_es, ref_es, method = "spearman", use = "pairwise.complete.obs"
        )),
        n_features = n(), .groups = "drop"
      ) %>%
      mutate(disease = disease_id, is_primary = gwas == primary_gwas)
  })
}

plot_data <- bind_rows(
  profile_consistency(spatial_profiles, "sr") %>% mutate(analysis = "spatial-region"),
  profile_consistency(cell_profiles, "cell_type") %>% mutate(analysis = "cell-type")
) %>%
  mutate(
    disease = factor(disease, levels = disease_order),
    analysis = factor(analysis, levels = c("spatial-region", "cell-type")),
    gwas_type = if_else(is_primary, "Primary GWAS", "Alternative GWAS")
  )
write_table(spatial_profiles, "human_dlpfc_all_gwas_spatial_profiles.tsv")
write_table(cell_profiles, "human_dlpfc_all_gwas_cell_type_profiles.tsv")
write_table(plot_data, "human_dlpfc_gwas_profile_consistency.tsv")

theme_plot <- theme_classic(base_size = 22, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    axis.title = element_text(size = 25, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    strip.text = element_text(size = 22, color = "black"),
    strip.background = element_rect(fill = "grey94", color = "grey82", linewidth = 0.25),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18, color = "black"),
    plot.margin = margin(3, 8, 3, 8)
  )

disease_colors <- c(
  ADHD = "#D89000", ANX = "#4B9FD3", BIP = "#00906E",
  MDD = "#C7A000", SCZ = "#C06A9F"
)
plot_data <- plot_data %>%
  mutate(
    point_shape = if_else(gwas_type == "Primary GWAS", "Primary", "Alternative"),
    analysis = recode(as.character(analysis),
                      "spatial-region" = "spatial-region",
                      "cell-type" = "cell-type"),
    analysis = factor(analysis, levels = c("spatial-region", "cell-type"))
  )

p <- ggplot(plot_data, aes(disease, spearman_rho, color = disease, shape = point_shape)) +
  geom_hline(yintercept = 0.7, linewidth = 0.35, linetype = "dashed", color = "#777777") +
  geom_point(size = 4.2, alpha = 0.88,
             position = position_jitter(width = 0.13, height = 0, seed = 7)) +
  facet_grid(
    analysis ~ .,
    scales = "fixed",
    labeller = as_labeller(c(
      "spatial-region" = "spatial\nregion",
      "cell-type"      = "cell\ntype"
    ))
  ) +
  scale_color_manual(values = disease_colors, guide = "none") +
  scale_shape_manual(
    values = c(Primary = 18, Alternative = 16),
    breaks = c("Primary", "Alternative"),
    labels = c("Primary GWAS", "Alternative GWAS"), name = NULL
  ) +
  scale_y_continuous(limits = c(0.65, 1.03), breaks = c(0.7, 0.8, 0.9, 1.0)) +
  labs(x = NULL, y = "Spearman rho") + theme_plot +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 21),
    strip.text.y.right = element_text(
      size = 16,
      angle = 0,
      lineheight = 0.95
    ),
    
    legend.position = "inside",
    legend.position.inside = c(0.36, 0.73),
    legend.justification = c(1, 1),
    
    legend.text = element_text(size = 16),
    legend.key.size = unit(0.65, "lines"),
    legend.background = element_rect(
      fill = scales::alpha("white", 0.75),
      color = NA
    ),
    legend.margin = margin(3, 4, 3, 4)
  )

stem <- file.path(out_dir, "human_dlpfc_gwas_profile_consistency")
ggsave(paste0(stem, ".png"), p, width = 8.5, height = 5.8, dpi = 600, bg = "white")
ggsave(paste0(stem, ".pdf"), p, width = 8.5, height = 5.8,
       device = grDevices::cairo_pdf, bg = "white")
