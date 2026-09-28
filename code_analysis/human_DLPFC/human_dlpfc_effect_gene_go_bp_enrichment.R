suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(jsonlite)
  library(ggplot2)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")
project_root <- normalizePath(getwd())
result_root <- file.path(project_root, "results")
out_dir <- file.path(
  project_root,
  "plot",
  "human_dlpfc_effect_gene_go_bp_enrichment"
)
proxy_url <- "http://127.0.0.1:7851"
curl_bin <- "/usr/bin/curl"
gprofiler_url <- "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


all_disease_order <- c("ADHD", "ANX", "BIP", "MDD", "SCZ")
display_disease_order <- c("ANX", "BIP", "MDD", "SCZ")
shared_group <- "Shared by >=2 diseases"
group_order <- c(display_disease_order, shared_group)
samples <- paste0("V", 1:12)

min_overlap_genes <- 5L
max_fdr <- 0.05
min_fold_enrichment <- 1


candidate_terms_per_group <- 20L
representative_terms_to_plot <- 15L
semantic_similarity_cutoff <- 0.70
fdr_colour_cap <- 10

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

thresholds <- read_tsv(
  file.path(project_root, "pvalues.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    disease,
    gwas,
    sample = v,
    threshold = as.numeric(p_value)
  )

gene_data <- imap_dfr(selected_gwas, function(gwas, disease) {
  map_dfr(samples, function(sample) {
    path <- file.path(
      result_root,
      disease,
      gwas,
      sample,
      sprintf(
        "%s_gene_marker_score.feather.genes.hg38.condi.assoc.tsv",
        sample
      )
    )
    x <- read_tsv(path, show_col_types = FALSE)
    names(x)[1] <- "gene"

    x %>%
      transmute(
        gene = as.character(gene),
        conditional_p = as.numeric(Condi.ECS.P)
      ) %>%
      mutate(
        disease = disease,
        gwas = gwas,
        sample = sample
      )
  })
}) %>%
  left_join(thresholds, by = c("disease", "gwas", "sample")) %>%
  filter(!is.na(gene), gene != ".", gene != "")

tested_genes <- gene_data %>%
  distinct(disease, gene)

effect_genes <- gene_data %>%
  filter(
    is.finite(conditional_p),
    is.finite(threshold),
    conditional_p < threshold
  ) %>%
  distinct(disease, gwas, gene, sample) %>%
  count(disease, gwas, gene, name = "n_sig_samples") %>%
  filter(n_sig_samples == length(samples)) %>%
  distinct(disease, gene)

# "Shared" means present in at least two diseases, not necessarily all five.
shared_background <- tested_genes %>%
  count(gene) %>%
  filter(n >= 2) %>%
  pull(gene)

shared_effect <- effect_genes %>%
  count(gene) %>%
  filter(n >= 2) %>%
  pull(gene)

gene_sets <- c(
  setNames(
    lapply(
      all_disease_order,
      function(x) filter(effect_genes, disease == x)$gene
    ),
    all_disease_order
  ),
  setNames(list(shared_effect), shared_group)
)

backgrounds <- c(
  setNames(
    lapply(
      all_disease_order,
      function(x) filter(tested_genes, disease == x)$gene
    ),
    all_disease_order
  ),
  setNames(list(shared_background), shared_group)
)

write_table(
  effect_genes,
  "human_dlpfc_go_input_genes_all12.tsv"
)

empty_enrichment <- function() {
  tibble(
    Group = character(),
    GO_ID = character(),
    GO_Term = character(),
    Background_Genes = integer(),
    Effect_Genes = integer(),
    Term_Genes = integer(),
    Overlap_Genes = integer(),
    Gene_Ratio = double(),
    Fold_Enrichment = double(),
    Odds_Ratio = double(),
    FDR = double(),
    Significant = logical(),
    Overlap_Gene_Symbols = character()
  )
}

gprofiler_bp <- function(group) {
  query <- unique(gene_sets[[group]])
  background <- unique(backgrounds[[group]])

  if (length(query) == 0L) {
    message("Skipping empty all-12-section gene set: ", group)
    return(empty_enrichment())
  }

  body <- jsonlite::toJSON(
    list(
      organism = "hsapiens",
      query = I(as.list(query)),
      sources = I(list("GO:BP")),
      user_threshold = 1,
      significance_threshold_method = "fdr",
      domain_scope = "custom",
      background = I(as.list(background)),
      no_evidences = FALSE,
      ordered = FALSE
    ),
    auto_unbox = TRUE
  )

  body_file <- tempfile(fileext = ".json")
  writeLines(body, body_file, useBytes = TRUE)

  args <- c(
    "-k",
    "--proxy", proxy_url,
    "--silent",
    "--show-error",
    "--retry", "3",
    "--connect-timeout", "60",
    "--max-time", "600",
    "-H", shQuote("Content-Type: application/json"),
    "--data-binary",
    shQuote(paste0("@", normalizePath(body_file, winslash = "/"))),
    gprofiler_url
  )

  response <- system2(
    curl_bin,
    args,
    stdout = TRUE,
    stderr = TRUE
  )
  unlink(body_file)

  parsed <- jsonlite::fromJSON(
    paste(response, collapse = ""),
    simplifyVector = FALSE
  )

  if (!is.null(parsed$message) && is.null(parsed$result)) {
    stop(parsed$message)
  }

  result <- parsed$result
  if (length(result) == 0L) {
    return(empty_enrichment())
  }

  gene_meta <- parsed$meta$genes_metadata$query$query_1
  ensg_to_symbol <- unlist(map(
    names(gene_meta$mapping),
    function(symbol) {
      setNames(
        rep(symbol, length(gene_meta$mapping[[symbol]])),
        unlist(gene_meta$mapping[[symbol]])
      )
    }
  ))
  mapped_ensgs <- unlist(gene_meta$ensgs)

  map_dfr(result, function(x) {
    overlap <- unique(unname(
      ensg_to_symbol[mapped_ensgs[lengths(x$intersections) > 0]]
    ))

    odds_ratio <- (
      (x$intersection_size + 0.5) *
        (
          x$effective_domain_size - x$term_size - x$query_size +
            x$intersection_size + 0.5
        )
    ) / (
      (x$query_size - x$intersection_size + 0.5) *
        (x$term_size - x$intersection_size + 0.5)
    )

    tibble(
      Group = group,
      GO_ID = x$native,
      GO_Term = x$name,
      Background_Genes = x$effective_domain_size,
      Effect_Genes = x$query_size,
      Term_Genes = x$term_size,
      Overlap_Genes = x$intersection_size,
      Gene_Ratio = x$intersection_size / x$query_size,
      Fold_Enrichment =
        (x$intersection_size / x$query_size) /
        (x$term_size / x$effective_domain_size),
      Odds_Ratio = odds_ratio,
      FDR = x$p_value,
      Significant = x$p_value < max_fdr,
      Overlap_Gene_Symbols = paste(overlap, collapse = ",")
    )
  })
}

# GO enrichment
enrichment <- map_dfr(group_order, gprofiler_bp) %>%
  arrange(
    factor(Group, levels = group_order),
    FDR,
    desc(Fold_Enrichment)
  )

write_table(
  enrichment,
  "human_dlpfc_effect_gene_go_bp_enrichment.tsv"
)

# Build a broad candidate pool
display_candidates <- enrichment %>%
  filter(
    Overlap_Genes >= min_overlap_genes,
    Fold_Enrichment > min_fold_enrichment,
    FDR < max_fdr
  ) %>%
  mutate(
    Neg_Log10_FDR = -log10(pmax(FDR, 1e-300)),
    Selection_Score = Neg_Log10_FDR * log2(Fold_Enrichment)
  )

if (nrow(display_candidates) == 0L) {
  stop("No GO BP terms passed the display criteria.")
}

# The per-group limit is only used to control the candidate pool before
# redundancy reduction; it is not the final set displayed for each disease.
candidate_pool <- display_candidates %>%
  group_by(Group) %>%
  slice_max(
    order_by = Selection_Score,
    n = candidate_terms_per_group,
    with_ties = FALSE
  ) %>%
  ungroup()

# Once a term enters the candidate pool, use all groups in which that term
# passed the display criteria when estimating its cross-disease support.
candidate_ids <- candidate_pool %>%
  distinct(GO_ID)

candidate_summary <- display_candidates %>%
  semi_join(candidate_ids, by = "GO_ID") %>%
  group_by(GO_ID, GO_Term) %>%
  summarise(
    Disease_Count = n_distinct(Group[Group %in% display_disease_order]),
    Shared_Present = any(Group == shared_group),
    Best_Score = max(Selection_Score, na.rm = TRUE),
    Best_FDR = min(FDR, na.rm = TRUE),
    Best_Fold_Enrichment = max(Fold_Enrichment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Give moderate preference to terms occurring in several diseases,
    # while retaining strongly enriched disease-specific biology.
    Priority_Score = Best_Score * (1 + 0.35 * Disease_Count)
  )


go_ids <- candidate_summary$GO_ID

if (length(go_ids) == 1L) {
  candidate_summary$Redundancy_Cluster <- 1L
} else {
  similarity_matrix <- rrvgo::calculateSimMatrix(
    go_ids,
    orgdb = "org.Hs.eg.db",
    ont = "BP",
    method = "Rel"
  )

  similarity_matrix[!is.finite(similarity_matrix)] <- 0
  similarity_matrix <- pmin(pmax(similarity_matrix, 0), 1)
  diag(similarity_matrix) <- 1

  redundancy_tree <- hclust(
    as.dist(1 - similarity_matrix),
    method = "average"
  )
  redundancy_cluster <- cutree(
    redundancy_tree,
    h = 1 - semantic_similarity_cutoff
  )

  candidate_summary$Redundancy_Cluster <- unname(
    redundancy_cluster[candidate_summary$GO_ID]
  )
}

representative_terms <- candidate_summary %>%
  arrange(
    desc(Disease_Count),
    desc(Shared_Present),
    desc(Priority_Score),
    Best_FDR
  ) %>%
  group_by(Redundancy_Cluster) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(
    desc(Disease_Count),
    desc(Shared_Present),
    desc(Priority_Score),
    Best_FDR
  ) %>%
  slice_head(n = representative_terms_to_plot) %>%
  mutate(
    Display_Order = row_number(),
    Term_Label = if_else(
      duplicated(GO_Term) | duplicated(GO_Term, fromLast = TRUE),
      paste0(GO_Term, " [", GO_ID, "]"),
      GO_Term
    ),
    Term_Label = str_wrap(Term_Label, width = 38)
  )

write_table(
  candidate_summary %>%
    left_join(
      representative_terms %>%
        select(GO_ID, Selected = Display_Order),
      by = "GO_ID"
    ) %>%
    arrange(Redundancy_Cluster, desc(Priority_Score)),
  "human_dlpfc_go_bp_semantic_reduction.tsv"
)

write_table(
  representative_terms,
  "human_dlpfc_go_bp_representative_terms.tsv"
)

term_levels <- representative_terms %>%
  arrange(desc(Display_Order)) %>%
  pull(Term_Label)

plot_data <- display_candidates %>%
  semi_join(
    representative_terms %>% select(GO_ID),
    by = "GO_ID"
  ) %>%
  left_join(
    representative_terms %>%
      select(GO_ID, Display_Order, Redundancy_Cluster, Term_Label),
    by = "GO_ID"
  ) %>%
  mutate(
    Group = factor(Group, levels = group_order),
    Term_Label = factor(Term_Label, levels = term_levels),
    Neg_Log10_FDR_Capped = pmin(Neg_Log10_FDR, fdr_colour_cap)
  )

write_table(
  plot_data %>%
    select(
      Group,
      GO_ID,
      GO_Term,
      Redundancy_Cluster,
      Background_Genes,
      Effect_Genes,
      Term_Genes,
      Overlap_Genes,
      Gene_Ratio,
      Fold_Enrichment,
      Odds_Ratio,
      FDR,
      Neg_Log10_FDR,
      Overlap_Gene_Symbols
    ),
  "human_dlpfc_go_bp_matrix_plot_data.tsv"
)

missing_groups <- setdiff(group_order, unique(as.character(plot_data$Group)))
if (length(missing_groups) > 0L) {
  message(
    "No GO terms passed the display criteria for: ",
    paste(missing_groups, collapse = ", "),
    ". Their matrix columns are retained and left blank."
  )
}

percent_labels <- function(x) {
  paste0(formatC(100 * x, format = "fg", digits = 2), "%")
}

group_axis_labels <- setNames(group_order, group_order)
group_axis_labels[shared_group] <- "Shared\n(>=2 diseases)"

p <- ggplot(
  plot_data,
  aes(x = Group, y = Term_Label)
) +
  geom_vline(
    xintercept = length(display_disease_order) + 0.5,
    colour = "#777777",
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      size = Gene_Ratio,
      colour = Neg_Log10_FDR_Capped
    ),
    alpha = 0.92
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = group_axis_labels
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_colour_gradient(
    low = "#F6D7CF",
    high = "#A71930",
    name = expression(-log[10](FDR)),
    limits = c(-log10(max_fdr), fdr_colour_cap),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(2.8, 9.5),
    name = "Gene ratio",
    labels = percent_labels
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 15,
    base_family = "Arial"
  ) +
  theme(
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 18,
      angle = 0,
      hjust = 0.5,
      vjust = 1
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 18
    ),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.margin = margin(10, 12, 10, 10)
  ) +
  guides(
    colour = guide_colourbar(
      order = 1,
      barheight = grid::unit(45, "mm"),
      barwidth = grid::unit(5, "mm")
    ),
    size = guide_legend(
      order = 2,
      override.aes = list(colour = "#A71930", alpha = 0.92)
    )
  )

stem <- file.path(
  out_dir,
  "human_dlpfc_effect_gene_go_bp_comparative_dotplot"
)

ggsave(
  paste0(stem, ".png"),
  p,
  width = 12.5,
  height = 9.5,
  dpi = 400,
  bg = "white"
)

grDevices::cairo_pdf(
  paste0(stem, ".pdf"),
  width = 12.5,
  height = 9.5,
  bg = "white"
)
print(p)
grDevices::dev.off()

message("Finished. Figure and tables were written to: ", out_dir)
