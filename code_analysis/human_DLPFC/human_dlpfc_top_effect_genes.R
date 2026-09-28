suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggrepel)
  library(grid)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")
project_root <- normalizePath(getwd())
result_root <- file.path(project_root, "results")
out_dir <- file.path(project_root, "plot", "human_dlpfc_top_effect_genes")
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

write_table <- function(x, name) {
  write_tsv(x, file.path(out_dir, name), na = "NA")
  invisible(x)
}

read_gene_assoc <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  names(x)[1] <- "gene"

  if ("#Var" %in% names(x)) {
    names(x)[names(x) == "#Var"] <- "n_var"
  }

  x %>%
    mutate(
      gene = as.character(gene),
      Condi.ECS.P = as.numeric(Condi.ECS.P),
      ECS.P = as.numeric(ECS.P),
      Score = if ("Score" %in% names(.)) {
        as.numeric(Score)
      } else {
        NA_real_
      }
    )
}

thresholds <- read_tsv(
  file.path(project_root, "pvalues.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    disease = as.character(disease),
    gwas = as.character(gwas),
    sample = as.character(v),
    condi_ecs_p_threshold = as.numeric(p_value)
  )

gene_data <- imap_dfr(selected_gwas, function(gwas_name, disease_id) {
  map_dfr(samples, function(sample_id) {
    path <- file.path(
      result_root,
      disease_id,
      gwas_name,
      sample_id,
      sprintf(
        "%s_gene_marker_score.feather.genes.hg38.condi.assoc.tsv",
        sample_id
      )
    )

    read_gene_assoc(path) %>%
      mutate(
        disease = disease_id,
        gwas = gwas_name,
        sample = sample_id
      )
  })
}) %>%
  left_join(sample_pairs, by = "sample") %>%
  left_join(
    thresholds,
    by = c("disease", "gwas", "sample")
  ) %>%
  mutate(
    is_effect_gene =
      is.finite(Condi.ECS.P) &
      is.finite(condi_ecs_p_threshold) &
      Condi.ECS.P < condi_ecs_p_threshold
  )

# A gene must pass the existing section-specific threshold in every V1-V12
# section. Missing/non-finite P values cannot count as significant, and
# duplicate rows cannot increase the number of significant sections.
stable_genes <- gene_data %>%
  filter(
    is_effect_gene,
    !is.na(gene),
    gene != ".",
    gene != ""
  ) %>%
  distinct(disease, gwas, gene, sample) %>%
  count(
    disease,
    gwas,
    gene,
    name = "n_sig_samples"
  ) %>%
  filter(n_sig_samples == length(samples))

effect_gene_data <- gene_data %>%
  filter(
    is_effect_gene,
    !is.na(gene),
    gene != ".",
    gene != ""
  ) %>%
  semi_join(
    stable_genes,
    by = c("disease", "gwas", "gene")
  ) %>%
  transmute(
    Disease = disease,
    GWAS = gwas,
    Tissue_Section = sample,
    Section_Pair = pair,
    Gene = gene,
    Ensembl_ID = EnsemblID,
    Chromosome,
    Start_Position = StartPosition,
    End_Position = EndPosition,
    N_Variants = n_var,
    ECS_P = ECS.P,
    Conditional_ECS_P = Condi.ECS.P,
    Conditional_ECS_P_Threshold = condi_ecs_p_threshold,
    Score
  ) %>%
  arrange(
    factor(Disease, levels = disease_order),
    Tissue_Section,
    Conditional_ECS_P,
    Gene
  )

write_table(
  effect_gene_data,
  "human_dlpfc_all_effect_genes.tsv"
)

effect_gene_summary <- effect_gene_data %>%
  group_by(
    Disease,
    GWAS,
    Gene,
    Ensembl_ID,
    Chromosome,
    Start_Position,
    End_Position
  ) %>%
  summarise(
    N_Significant_Sections = n_distinct(Tissue_Section),
    Significant_Sections = paste(
      sort(unique(Tissue_Section)),
      collapse = ","
    ),
    Min_Conditional_ECS_P = min(Conditional_ECS_P),
    Median_Conditional_ECS_P = median(Conditional_ECS_P),
    Max_Score = if (all(is.na(Score))) {
      NA_real_
    } else {
      max(Score, na.rm = TRUE)
    },
    Median_Score = if (all(is.na(Score))) {
      NA_real_
    } else {
      median(Score, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  arrange(
    factor(Disease, levels = disease_order),
    desc(N_Significant_Sections),
    Min_Conditional_ECS_P,
    Gene
  )

write_table(
  effect_gene_summary,
  "human_dlpfc_effect_gene_summary.tsv"
)

rank_rows <- gene_data %>%
  filter(
    is_effect_gene,
    !is.na(gene),
    gene != ".",
    gene != ""
  ) %>%
  semi_join(
    stable_genes,
    by = c("disease", "gwas", "gene")
  )

loc_genes_excluded_from_plot <- rank_rows %>%
  filter(grepl("^LOC[0-9]+$", gene)) %>%
  distinct(disease, gwas, gene) %>%
  arrange(factor(disease, levels = disease_order), gene)

write_table(
  loc_genes_excluded_from_plot,
  "human_dlpfc_loc_genes_excluded_from_plot.tsv"
)

plot_data <- rank_rows %>%
  semi_join(
    stable_genes,
    by = c("disease", "gene")
  ) %>%
  # Exclude provisional LOC-plus-digits symbols from the displayed ranking.
  # They remain present in the complete effect-gene result tables above.
  filter(!grepl("^LOC[0-9]+$", gene)) %>%
  group_by(disease, sample) %>%
  arrange(Condi.ECS.P, .by_group = TRUE) %>%
  mutate(sample_rank = row_number()) %>%
  ungroup() %>%
  group_by(disease, gene) %>%
  summarise(
    mean_rank = mean(sample_rank, na.rm = TRUE),
    min_condi_p = min(Condi.ECS.P, na.rm = TRUE),
    median_condi_p = median(Condi.ECS.P, na.rm = TRUE),
    max_score = if (all(!is.finite(Score))) {
      NA_real_
    } else {
      max(Score[is.finite(Score)])
    },
    median_score = median(Score, na.rm = TRUE),
    mean_score = mean(Score, na.rm = TRUE),
    n_sig_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  group_by(disease) %>%
  arrange(
    mean_rank,
    min_condi_p,
    desc(max_score),
    .by_group = TRUE
  ) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  mutate(
    disease = factor(disease, levels = disease_order),
    neglog10_median_p = -log10(pmax(median_condi_p, 1e-300))
  ) %>%
  group_by(disease) %>%
  arrange(neglog10_median_p, .by_group = TRUE) %>%
  mutate(y_pos = row_number()) %>%
  ungroup()

# Shared genes are defined using the complete stable-effect-gene sets, not
# only the displayed Top 20 genes.
shared_gene_counts <- stable_genes %>%
  distinct(disease, gene) %>%
  count(gene, name = "n_diseases")

plot_data <- plot_data %>%
  left_join(shared_gene_counts, by = "gene") %>%
  mutate(
    Shared_Status = if_else(
      n_diseases >= 2,
      "Shared gene",
      "Disease-specific"
    ),
    Shared_Status = factor(
      Shared_Status,
      levels = c("Disease-specific", "Shared gene")
    )
  )

write_table(
  plot_data,
  "human_dlpfc_top_effect_genes.tsv"
)

disease_colors <- c(
  ADHD = "#D89000",
  ANX = "#4B9FD3",
  BIP = "#00906E",
  MDD = "#C7A000",
  SCZ = "#C06A9F"
)

# Very light versions of the disease colours, used only as panel backgrounds.
disease_background_colors <- c(
  ADHD = "#FFF8E8",
  ANX = "#EEF7FC",
  BIP = "#ECF8F4",
  MDD = "#FFFBEA",
  SCZ = "#FAF0F6"
)

panel_background <- tibble(
  disease = factor(disease_order, levels = disease_order)
)

theme_plot <- theme_classic(
  base_size = 22,
  base_family = "Arial"
) +
  theme(
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    axis.title = element_text(size = 42, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    strip.text = element_text(size = 32, face = "bold", color = "black"),
    strip.background = element_rect(
      fill = "grey94",
      color = "grey82",
      linewidth = 0.25
    ),
    panel.grid = element_blank(),
    plot.margin = margin(3, 3, 3, 3)
  )

p <- ggplot(
  plot_data,
  aes(neglog10_median_p, y_pos, color = disease)
) +
  geom_rect(
    data = panel_background,
    aes(fill = disease),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    inherit.aes = FALSE,
    color = NA
  ) +
  geom_hline(
    yintercept = seq(1, 20, 1),
    color = "grey88",
    linewidth = 0.18
  ) +
  geom_point(
    aes(
      size = ifelse(is.finite(median_score), median_score, 0),
      shape = Shared_Status
    ),
    alpha = 0.94
  ) +
  geom_text_repel(
    aes(label = gene),
    color = "black",
    fontface = "italic",
    size = 10,
    min.segment.length = 0,
    segment.size = 0.16,
    segment.alpha = 0.42,
    box.padding = 0.22,
    point.padding = 0.12,
    max.overlaps = Inf,
    force = 10,
    force_pull = 0.03,
    direction = "both",
    seed = 7
  ) +
  facet_wrap(
    ~disease,
    ncol = 1,
    scales = "free_x"
  ) +
  scale_color_manual(
    values = disease_colors,
    guide = "none"
  ) +
  scale_fill_manual(
    values = disease_background_colors,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "Disease-specific" = 16,
      "Shared gene" = 18
    ),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(4, 9),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.06, 0.18))
  ) +
  scale_y_continuous(
    limits = c(0.45, 20.60),
    breaks = NULL
  ) +
  labs(
    x = "−log10(median conditional ECS P)",
    y = "effect genes"
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        color = "#555555",
        size = 9,
        alpha = 1
      )
    )
  ) +
  theme_plot +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15, l = 15)),
    panel.spacing.y = unit(5, "pt"),
    strip.text = element_text(size = 22, hjust = 0),
    legend.position = c(0.985, 0.965),
    legend.justification = c(1, 1),
    legend.direction = "horizontal",
    legend.text = element_text(size = 36, color = "black"),
    legend.key.width = unit(0.9, "cm"),
    legend.key.height = unit(0.9, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    legend.background = element_rect(
      fill = "#FFFFFFE6",
      color = "grey65",
      linewidth = 0.4
    ),
    legend.margin = margin(5, 8, 5, 8)
  )

if (nrow(plot_data) == 0L) {
  p <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No genes significant in all 12 sections"
    ) +
    theme_void()
}

stem <- file.path(out_dir, "human_dlpfc_top_effect_genes")

ggsave(
  paste0(stem, ".png"),
  p,
  width = 40,
  height = 22,
  dpi = 300,
  bg = "white"
)

ggsave(
  paste0(stem, ".pdf"),
  p,
  width = 40,
  height = 22,
  device = grDevices::cairo_pdf,
  bg = "white"
)
