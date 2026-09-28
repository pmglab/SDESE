suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(grid)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")

project_root <- normalizePath(getwd())
result_root <- file.path(project_root, "results")
out_dir <- file.path(
  project_root,
  "plot",
  "human_dlpfc_driver_gene_jaccard"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

disease_order <- c("ADHD", "ANX", "BIP", "MDD", "SCZ")
samples <- paste0("V", 1:12)

selected_gwas <- c(
  ADHD = "ADHD2022_iPSYCH_deCODE_PGC",
  ANX = "ANX_2026_daner_fullANX_v12_woUTAH_11022026",
  BIP = "bip2024_eur_noUKB_no23andMe",
  MDD = "pgc_mdd2025_no23andMe_noUKBB_eur_v3_49_24_11",
  SCZ = "PGC3_SCZ_wave3_european_autosome_public_v3_vcf"
)

thresholds <- read_tsv(
  file.path(project_root, "pvalues.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    disease = as.character(disease),
    gwas = as.character(gwas),
    sample = as.character(v),
    threshold = as.numeric(p_value)
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

    x <- read_tsv(path, show_col_types = FALSE)
    names(x)[1] <- "gene"

    x %>%
      transmute(
        disease = disease_id,
        gwas = gwas_name,
        sample = sample_id,
        gene = as.character(gene),
        Score = as.numeric(Score),
        Condi.ECS.P = as.numeric(Condi.ECS.P)
      )
  })
}) %>%
  filter(
    !is.na(gene),
    gene != "",
    gene != "."
  ) %>%
  left_join(
    thresholds,
    by = c("disease", "gwas", "sample")
  )

all_genes <- gene_data %>%
  distinct(disease, sample, gene) %>%
  count(disease, gene, name = "n_samples") %>%
  filter(n_samples == length(samples)) %>%
  select(disease, gene)

conditional_genes <- gene_data %>%
  filter(
    is.finite(Score),
    is.finite(Condi.ECS.P),
    is.finite(threshold),
    Condi.ECS.P < threshold
  ) %>%
  distinct(disease, sample, gene) %>%
  count(disease, gene, name = "n_sig_samples") %>%
  filter(n_sig_samples == length(samples)) %>%
  select(disease, gene)

make_gene_sets <- function(data) {
  setNames(
    lapply(disease_order, function(disease_id) {
      data %>%
        filter(disease == disease_id) %>%
        pull(gene) %>%
        unique()
    }),
    disease_order
  )
}

calculate_jaccard <- function(set_a, set_b) {
  union_size <- length(union(set_a, set_b))

  if (union_size == 0) {
    return(0)
  }

  length(intersect(set_a, set_b)) / union_size
}

all_gene_sets <- make_gene_sets(all_genes)
conditional_sets <- make_gene_sets(conditional_genes)

disease_pairs <- combn(
  disease_order,
  2,
  simplify = FALSE
)

comparison_data <- map_dfr(disease_pairs, function(pair) {
  disease_a <- pair[1]
  disease_b <- pair[2]

  data.frame(
    pair = paste(disease_a, "vs", disease_b),
    all_genes = calculate_jaccard(
      all_gene_sets[[disease_a]],
      all_gene_sets[[disease_b]]
    ),
    conditional = calculate_jaccard(
      conditional_sets[[disease_a]],
      conditional_sets[[disease_b]]
    )
  )
}) %>%
  mutate(
    reduction = all_genes - conditional
  )

pair_order <- comparison_data$pair

plot_data <- bind_rows(
  comparison_data %>%
    transmute(
      pair,
      analysis = "All genes",
      jaccard = all_genes
    ),
  comparison_data %>%
    transmute(
      pair,
      analysis = "Conditional significant genes",
      jaccard = conditional
    )
) %>%
  mutate(
    pair = factor(pair, levels = pair_order),
    analysis = factor(
      analysis,
      levels = c(
        "All genes",
        "Conditional significant genes"
      )
    )
  )

max_jaccard <- max(plot_data$jaccard, na.rm = TRUE)

y_upper <- 0.26

mean_all <- mean(comparison_data$all_genes, na.rm = TRUE)
mean_conditional <- mean(comparison_data$conditional, na.rm = TRUE)

dodge_width <- 0.82

legend_inside_theme <- if (packageVersion("ggplot2") >= "3.5.0") {
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98)
  )
} else {
  theme(
    legend.position = c(0.02, 0.98)
  )
}

p <- ggplot(
  plot_data,
  aes(
    x = pair,
    y = jaccard,
    fill = analysis
  )
) +
  geom_col(
    position = position_dodge(width = dodge_width),
    width = 0.72,
    color = "white",
    linewidth = 0.45
  ) +
  geom_text(
    aes(label = sprintf("%.3f", jaccard)),
    position = position_dodge(width = dodge_width),
    vjust = -0.45,
    size = 4,
    family = "Arial",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "All genes" = "#377EB8",
      "Conditional significant genes" = "#E41A1C"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    breaks = pretty(c(0, y_upper), n = 6),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(
    ylim = c(0, y_upper),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Jaccard index",
    subtitle = sprintf(
      "Mean Jaccard: %.3f → %.3f",
      mean_all,
      mean_conditional
    )
  ) +
  theme_classic(
    base_size = 21,
    base_family = "Arial"
  ) +
  theme(
    axis.title.y = element_text(
      size = 25,
      face = "bold",
      color = "black",
      margin = margin(r = 12)
    ),
    axis.text.x = element_text(
      size = 18,
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.7
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 17, color = "black"),
    legend.key.height = unit(0.65, "cm"),
    legend.key.width = unit(0.9, "cm"),
    legend.background = element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = 0.88),
      color = "#666666",
      linewidth = 0.45
    ),
    legend.margin = margin(6, 8, 6, 8),
    plot.subtitle = element_text(
      size = 19,
      color = "black",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.margin = margin(12, 14, 10, 10)
  ) +
  legend_inside_theme

write_tsv(
  all_genes,
  file.path(out_dir, "human_dlpfc_all_recurrent_genes.tsv")
)

write_tsv(
  conditional_genes,
  file.path(out_dir, "human_dlpfc_conditional_recurrent_genes.tsv")
)

write_tsv(
  comparison_data,
  file.path(out_dir, "human_dlpfc_jaccard_comparison.tsv")
)

ggsave(
  file.path(out_dir, "human_dlpfc_jaccard_comparison_bar.png"),
  p,
  width = 10,
  height = 7,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(out_dir, "human_dlpfc_jaccard_comparison_bar.pdf"),
  p,
  width = 10,
  height = 7,
  device = grDevices::cairo_pdf,
  bg = "white",
  limitsize = FALSE
)
