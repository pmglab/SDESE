suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(readr)
  library(hdf5r)
  library(Matrix)
})

project_root <- normalizePath(getwd(), winslash = "/")

xlsx_file <- file.path(project_root, "41593_2020_787_MOESM3_ESM.xlsx")
h5ad_dir <- file.path(project_root, "resources", "Brain_33558695_h5ad")
out_dir <- file.path(project_root, "plot", "human_dlpfc_sr_layer_marker_profiles")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, file, ...) {
  readr::write_tsv(x, file, na = "NA", ...)
  invisible(x)
}

sr_order <- c("SR7", "SR5", "SR6", "SR1", "SR4", "SR2", "SR3")
layer_order <- c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM")
layer_labels <- c(
  Layer1 = "L1",
  Layer2 = "L2",
  Layer3 = "L3",
  Layer4 = "L4",
  Layer5 = "L5",
  Layer6 = "L6",
  WM = "WM"
)
sr_labels <- c(
  SR7 = "SR7",
  SR5 = "SR5",
  SR6 = "SR6",
  SR1 = "SR1",
  SR4 = "SR4",
  SR2 = "SR2",
  SR3 = "SR3"
)
layer_group_map <- c(
  Layer1 = "Superficial",
  Layer2 = "Superficial",
  Layer3 = "Middle",
  Layer4 = "Middle",
  Layer5 = "Deep",
  Layer6 = "Deep",
  WM = "WM"
)
sr_group_map <- c(
  SR7 = "Superficial",
  SR5 = "Superficial",
  SR6 = "Middle",
  SR1 = "Middle",
  SR4 = "Deep",
  SR2 = "Deep",
  SR3 = "WM"
)
laminar_group_order <- c("Superficial", "Middle", "Deep", "WM")

read_obs_column <- function(h5, column) {
  obs <- h5[["obs"]]
  x <- obs[[column]]
  if (inherits(x, "H5Group") && "codes" %in% names(x) && "categories" %in% names(x)) {
    codes <- as.integer(x[["codes"]]$read())
    categories <- as.character(x[["categories"]]$read())
    return(categories[codes + 1L])
  }
  as.character(x$read())
}

read_var_gene_names <- function(h5) {
  as.character(h5[["var"]][["gene"]]$read())
}

read_x_matrix <- function(h5) {
  x <- h5[["X"]]
  shape <- as.integer(x$attr_open("shape")$read())
  t(new(
    "dgCMatrix",
    x = as.numeric(x[["data"]]$read()),
    i = as.integer(x[["indices"]]$read()),
    p = as.integer(x[["indptr"]]$read()),
    Dim = as.integer(c(shape[2], shape[1]))
  ))
}

read_sample_sr_expression <- function(h5ad_file) {
  h5 <- H5File$new(h5ad_file, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  sample_id <- sub("^HS_Brain_33558695_(V[0-9]+)\\.seuratobj\\.h5ad$", "\\1", basename(h5ad_file))
  genes <- make.unique(read_var_gene_names(h5))
  sr <- paste0("SR", read_obs_column(h5, "spatial_region"))
  keep_cells <- sr %in% sr_order

  x <- read_x_matrix(h5)
  colnames(x) <- genes

  x <- x[keep_cells, , drop = FALSE]
  sr <- sr[keep_cells]
  libsize <- Matrix::rowSums(x)
  libsize[libsize == 0] <- 1
  x_norm <- Diagonal(x = 1e4 / libsize) %*% x
  x_norm@x <- log1p(x_norm@x)

  rows <- lapply(sr_order, function(sr_id) {
    cell_idx <- which(sr == sr_id)
    if (length(cell_idx) == 0) return(NULL)
    tibble(
      sample = sample_id,
      sr = sr_id,
      gene = colnames(x_norm),
      mean_expr = as.numeric(Matrix::colMeans(x_norm[cell_idx, , drop = FALSE])),
      n_spots = length(cell_idx)
    )
  })
  bind_rows(rows)
}

marker_top_n <- 100
high_z_threshold <- 0.9

marker_table <- read_excel(xlsx_file, sheet = "Table S4B") %>%
  mutate(gene = toupper(gene), ensembl = as.character(ensembl))
marker_long <- lapply(layer_order, function(layer_id) {
  marker_table %>%
    transmute(
      layer = layer_id, gene, ensembl,
      fdr = .data[[paste0("fdr_", layer_id)]],
      t_stat = .data[[paste0("t_stat_", layer_id)]]
    ) %>%
    filter(!is.na(gene), gene != "", is.finite(t_stat), t_stat > 0, is.finite(fdr)) %>%
    arrange(fdr, desc(t_stat)) %>%
    slice_head(n = marker_top_n)
}) %>%
  bind_rows() %>%
  distinct(layer, gene, .keep_all = TRUE)

h5ad_files <- list.files(h5ad_dir, pattern = "\\.seuratobj\\.h5ad$", full.names = TRUE)
sr_expr <- bind_rows(lapply(h5ad_files, read_sample_sr_expression)) %>%
  mutate(gene = toupper(gene)) %>%
  group_by(sample, sr, gene) %>%
  summarise(
    mean_expr = mean(mean_expr, na.rm = TRUE),
    n_spots = max(n_spots, na.rm = TRUE),
    .groups = "drop"
  )

sample_layer_scores <- sr_expr %>%
  inner_join(marker_long %>% select(layer, gene), by = "gene", relationship = "many-to-many") %>%
  group_by(sample, sr, layer) %>%
  summarise(
    marker_score = mean(mean_expr, na.rm = TRUE),
    n_marker_detected = n_distinct(gene),
    .groups = "drop"
  ) %>%
  group_by(sample, layer) %>%
  mutate(marker_score_z = as.numeric(scale(marker_score))) %>%
  ungroup()

mapping_scores <- sample_layer_scores %>%
  group_by(sr, layer) %>%
  summarise(
    mean_z = mean(marker_score_z, na.rm = TRUE),
    se_z = sd(marker_score_z, na.rm = TRUE) / sqrt(n_distinct(sample)),
    mean_detected_markers = mean(n_marker_detected, na.rm = TRUE),
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  mutate(
    sr = factor(sr, levels = sr_order),
    layer = factor(layer, levels = layer_order)
  )

write_table(marker_long, file.path(out_dir, "maynard_dlpfc_layer_marker_genes_top100.tsv"))
write_table(sample_layer_scores, file.path(out_dir, "human_dlpfc_sample_sr_layer_marker_scores.tsv"))
write_table(mapping_scores, file.path(out_dir, "human_dlpfc_sr_layer_marker_mapping_scores.tsv"))

profile_scores <- mapping_scores %>%
  mutate(
    sr_label = factor(sr_labels[as.character(sr)], levels = sr_labels[sr_order]),
    sr_index = match(as.character(sr), sr_order),
    layer_label = factor(layer_labels[as.character(layer)], levels = layer_labels[layer_order]),
    layer_group = factor(layer_group_map[as.character(layer)], levels = laminar_group_order),
    sr_group = factor(sr_group_map[as.character(sr)], levels = laminar_group_order)
  ) %>%
  group_by(layer_label) %>%
  mutate(
    is_highlight = mean_z >= high_z_threshold | mean_z == max(mean_z, na.rm = TRUE),
    text_color   = ifelse(abs(mean_z) > 1.15, "white", "black")
  ) %>%
  ungroup()


z_lim <- max(abs(profile_scores$mean_z), na.rm = TRUE) * 1.08

sr_group_label_pos <- tibble(
  group = factor(laminar_group_order, levels = laminar_group_order),
  x = c(1.5, 3.5, 5.5, 7),
  y = 7.62
)

gg <- ggplot(profile_scores, aes(x = sr_label, y = layer_label)) +

  geom_tile(aes(fill = mean_z), color = "white", linewidth = 0.8) +

  geom_tile(
    data = profile_scores %>% filter(is_highlight),
    fill = NA, color = "black", linewidth = 0.9
  ) +

  geom_text(
    aes(
      label = sprintf("%.2f", mean_z),
      fontface = ifelse(is_highlight, "bold", "plain"),
      color = text_color
    ),
    size = 3.2
  ) +

  geom_vline(xintercept = c(2.5, 4.5, 6.5), color = "grey45", linewidth = 0.4) +
  geom_hline(yintercept = c(2.5, 4.5, 6.5), color = "grey45", linewidth = 0.4) +

  annotate("text",
           x = sr_group_label_pos$x, y = sr_group_label_pos$y,
           label = as.character(sr_group_label_pos$group),
           size = 3, fontface = "bold", color = "grey30"
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0,
    breaks = seq(-2, 2, by = 1),
    limits = c(-z_lim, z_lim),
    name = "Marker expression\nz-score"
  ) +
  scale_color_identity() +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.6))) +
  scale_y_discrete(
    limits = rev(levels(profile_scores$layer_label)),
    expand = expansion(add = c(0.5, 1.25))
  ) +
  labs(x = "PSC spatial region", y = "Cortical layer") +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 10.5, face = "bold", margin = margin(t = 4)),
    axis.text.y = element_text(size = 10.5, face = "bold", margin = margin(r = 4)),
    axis.title.x = element_text(margin = margin(t = 8), size = 10.5),
    axis.title.y = element_text(margin = margin(r = 8), size = 10.5),
    legend.position = "right",
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.45, "in"),
    plot.margin = margin(10, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "human_dlpfc_sr_layer_marker_profiles.png"),
  gg, width = 6.6, height = 6.2, dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "human_dlpfc_sr_layer_marker_profiles.pdf"),
  gg, width = 6.6, height = 6.2, device = grDevices::cairo_pdf, bg = "white"
)
