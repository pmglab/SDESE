suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

setwd("/public3/ly/SDESE_summary/human_DLPFC")
project_dir <- normalizePath(getwd())
results_dir <- file.path(project_dir, "results")
pvalues_file <- file.path(project_dir, "pvalues.tsv")
out_dir <- Sys.getenv("SDESE_PLOT_OUTPUT", unset = file.path(project_dir, "plot", "sdese_effect_gene_gwas_robustness"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- paste0("V", 1:12)
fdr_threshold <- 0.05
gwas_labels <- c(
  ADHD2022_iPSYCH_deCODE_PGC = "ADHD-2022",
  ADHD_OVERALL_vanderlaan2025_GWAS_summary_stats = "ADHD-2025",
  ANX_2026_daner_fullANX_v12_woUTAH_11022026 = "ANX-2026",
  pgc_anx2_gadsymptsquant_2026_eur = "PGC-ANX2",
  bip2024_eur_noUKB_no23andMe = "BIP-2024",
  daner_bip_pgc3_nm_noukbiobank = "PGC3-BIP",
  pgc_bip2021_all_vcf = "BIP-2021",
  pgc_mdd2025_no23andMe_noUKBB_eur_v3_49_24_11 = "MDD-2025",
  `PGC_UKB_depression_genome-wide` = "PGC-UKB-MDD",
  CLOZUK_PGC2noclo_METAL_assoc_dosage = "CLOZUK-PGC2",
  PGC3_SCZ_wave3_european_autosome_public_v3_vcf = "PGC3-SCZ",
  daner_PGC_SCZ52_0513a_hq2 = "PGC2-SCZ52",
  daner_natgen_pgc_eur = "SCZ-NatGen"
)
disease_colors <- c(
  ADHD = "#3775BA", ANX = "#42949E", BIP = "#9A4D8E",
  MDD = "#B86B77", SCZ = "#5C647E"
)

# Use the pipeline's existing section-specific BH-derived cutoffs.
threshold_table <- read_tsv(pvalues_file, show_col_types = FALSE) %>%
  transmute(disease = as.character(disease), gwas = as.character(gwas),
            sample = as.character(v), threshold = as.numeric(p_value))
if (anyDuplicated(threshold_table[c("disease", "gwas", "sample")])) {
  stop("Duplicate disease/GWAS/section entries in pvalues.tsv")
}

read_gene_file <- function(path) {
  if (!file.exists(path)) stop("Missing required section file: ", path)
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(gene = as.character(RegionID), condi_p = as.numeric(Condi.ECS.P)) %>%
    filter(!is.na(gene), gene != "", gene != ".") %>%
    group_by(gene) %>%
    summarise(condi_p = if (any(is.finite(condi_p))) min(condi_p[is.finite(condi_p)]) else NA_real_,
              .groups = "drop")
}

load_gwas <- function(gwas_dir, disease) {
  gwas <- basename(gwas_dir)
  rows <- map_dfr(samples, function(sample) {
    cutoff <- threshold_table %>%
      filter(.data$disease == .env$disease, .data$gwas == .env$gwas,
             .data$sample == .env$sample) %>% pull(threshold)
    if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff > 1) {
      stop("Missing or invalid threshold: ", disease, "/", gwas, "/", sample)
    }
    path <- file.path(gwas_dir, sample,
      sprintf("%s_gene_marker_score.feather.genes.hg38.condi.assoc.tsv", sample))
    read_gene_file(path) %>% mutate(sample = sample, threshold = cutoff,
      significant = is.finite(condi_p) & condi_p < cutoff)
  })
  # Non-finite/missing results never count toward all-12 significance.
  summary <- rows %>% group_by(gene) %>% summarise(
    N_Tested_Sections = n_distinct(sample[is.finite(condi_p)]),
    N_Significant_Sections = n_distinct(sample[significant]), .groups = "drop")
  list(tested = summary$gene[summary$N_Tested_Sections == length(samples)],
       significant = summary$gene[summary$N_Significant_Sections == length(samples)],
       summary = summary %>% mutate(Disease = disease, GWAS = gwas))
}

compare_sets <- function(data1, data2) {
  # Pair-specific universe: genes with valid conditional tests in all 12
  # sections in BOTH GWAS. No additional marginal-ECS or Score filter.
  background <- intersect(data1$tested, data2$tested)
  set1 <- intersect(data1$significant, background)
  set2 <- intersect(data2$significant, background)
  population <- length(background)
  n1 <- length(set1); n2 <- length(set2)
  overlap <- length(intersect(set1, set2))
  expected <- if (population > 0L) n1 * n2 / population else NA_real_
  tibble(
    N_Required_Sections = length(samples),
    Background_Common_Tested_All12_Genes = population,
    GWAS1_All12_Significant_Genes = length(data1$significant),
    GWAS2_All12_Significant_Genes = length(data2$significant),
    GWAS1_Significant_In_Background = n1,
    GWAS2_Significant_In_Background = n2,
    Observed_Overlap = overlap, Expected_Overlap = expected,
    Fold_Enrichment = if (is.finite(expected) && expected > 0) overlap / expected else NA_real_,
    Hypergeometric_P = if (population == 0L || n1 == 0L || n2 == 0L) NA_real_ else
      phyper(overlap - 1L, n1, population - n1, n2, lower.tail = FALSE),
    Overlap_Genes = paste(sort(intersect(set1, set2)), collapse = ",")
  )
}

all_gene_summaries <- list()
analyse_disease <- function(disease_dir) {
  disease <- basename(disease_dir)
  gwas_dirs <- list.dirs(disease_dir, recursive = FALSE, full.names = TRUE)
  if (length(gwas_dirs) < 2L) stop("Fewer than two GWAS datasets for ", disease)
  gwas_data <- setNames(map(gwas_dirs, load_gwas, disease = disease), basename(gwas_dirs))
  all_gene_summaries[[disease]] <<- map_dfr(gwas_data, "summary")
  map_dfr(combn(names(gwas_data), 2, simplify = FALSE), function(pair) {
    compare_sets(gwas_data[[pair[1]]], gwas_data[[pair[2]]]) %>%
      mutate(Disease = disease, GWAS1 = pair[1], GWAS2 = pair[2], .before = 1)
  })
}

detail <- map_dfr(file.path(results_dir, names(disease_colors)), analyse_disease) %>%
  mutate(
    # One overlap test per GWAS pair; BH correction across all plotted pairs.
    Hypergeometric_FDR = p.adjust(Hypergeometric_P, method = "BH"),
    Hypergeometric_FDR_Significant = case_when(
      is.na(Hypergeometric_FDR) ~ "Not testable",
      Hypergeometric_FDR < fdr_threshold ~ "Yes", TRUE ~ "No")
  ) %>% arrange(Disease, GWAS1, GWAS2)

summary_data <- detail %>% mutate(
  GWAS_Pair = paste(coalesce(unname(gwas_labels[GWAS1]), GWAS1), "vs",
                    coalesce(unname(gwas_labels[GWAS2]), GWAS2)),
  FDR_Label = if_else(is.na(Hypergeometric_FDR), "NA",
                     formatC(Hypergeometric_FDR, format = "e", digits = 1)),
  y = rev(seq_len(n()))
)
gene_summary <- bind_rows(all_gene_summaries)
write_tsv(gene_summary %>% filter(N_Significant_Sections == length(samples)),
          file.path(out_dir, "all_gwas_effect_genes_all12.tsv"))
write_tsv(detail, file.path(out_dir, "all_diseases_gwas_hypergeometric.tsv"), na = "NA")
write_tsv(summary_data, file.path(out_dir, "sdese_effect_gene_gwas_robustness_summary.tsv"), na = "NA")

# One point is one all-12-section gene-set comparison, not 12 replicates.
x_max <- max(c(1, summary_data$Fold_Enrichment), na.rm = TRUE)
annotation_x <- x_max * 1.22
boundaries <- summary_data %>% mutate(next_disease = lead(Disease), next_y = lead(y)) %>%
  filter(Disease != next_disease) %>% transmute(y = (y + next_y) / 2)
p <- ggplot(summary_data, aes(Fold_Enrichment, y)) +
  geom_hline(data = boundaries, aes(yintercept = y), inherit.aes = FALSE,
             color = "#D9D9D9", linewidth = 0.4) +
  geom_vline(xintercept = 1, color = "#767676", linewidth = 0.7, linetype = "dashed") +
  geom_point(aes(fill = Disease, shape = Hypergeometric_FDR_Significant),
             color = "#444444", size = 6, stroke = 0.6, na.rm = TRUE) +
  geom_text(aes(x = annotation_x, label = FDR_Label), hjust = 0.5,
            size = 8.5, color = "#272727") +
  annotate("text", x = annotation_x, y = nrow(summary_data) + 0.85,
           label = "Overlap test\nBH-FDR", size = 8.5, lineheight = 0.95) +
  scale_fill_manual(values = disease_colors, guide = "none") +
  scale_shape_manual(values = c("Yes" = 21, "No" = 1, "Not testable" = 4), guide = "none") +
  scale_y_continuous(breaks = summary_data$y, labels = summary_data$GWAS_Pair,
                     limits = c(0.5, nrow(summary_data) + 1.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, x_max * 1.45),
                     breaks = pretty(c(0, x_max), n = 5), expand = c(0, 0)) +
  labs(x = "Fold enrichment of overlapping effect genes", y = NULL) +
  theme_classic(base_size = 22, base_family = "Arial") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 26, margin = margin(r = 12)),
        axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 32),
        plot.margin = margin(12, 24, 12, 12))
height <- max(9.5, 0.68 * nrow(summary_data) + 3)
stem <- file.path(out_dir, "sdese_effect_gene_gwas_robustness")
ggsave(paste0(stem, ".png"), p, width = 15.5, height = height, dpi = 600, bg = "white")
ggsave(paste0(stem, ".pdf"), p, width = 15.5, height = height,
       device = grDevices::cairo_pdf, bg = "white")
message("All-12-section comparisons: ", nrow(detail), "; BH-FDR < 0.05: ",
        sum(detail$Hypergeometric_FDR < 0.05, na.rm = TRUE))
