suppressPackageStartupMessages({
  library(hdf5r)
  library(ggplot2)
  library(data.table)
  library(scales)
  library(ragg)
})

options(stringsAsFactors = FALSE)

input_h5ad <- "E16.5_E1S1.MOSTA.h5ad"
output_dir <- "heatmaps_fdr_bh"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

traits <- c("ADHD", "SCZ", "BIP", "MDD", "ANX", "IQ", "LC", "TP", "DBP", "CAD")

# Kept identical to the tissue order in the supplied reference heatmap.
tissues <- c(
  "Muscle", "Liver", "Brain", "Adrenal gland", "Inner ear", "Smooth muscle",
  "Lung", "Submandibular gland", "Spinal cord", "Cartilage primordium",
  "Connective tissue", "Sympathetic nerve", "Kidney", "Dorsal root ganglion",
  "GI tract", "Cartilage", "Mucosal epithelium", "Choroid plexus",
  "Adipose tissue", "Jaw and tooth", "Cavity", "Bone", "Heart", "Meninges",
  "Epidermis"
)

tissue_group <- c(
  "Muscle" = "Musculoskeletal", "Liver" = "Visceral organ", "Brain" = "Neural",
  "Adrenal gland" = "Visceral organ", "Inner ear" = "Neural",
  "Smooth muscle" = "Musculoskeletal", "Lung" = "Cardiopulmonary",
  "Submandibular gland" = "Visceral organ", "Spinal cord" = "Neural",
  "Cartilage primordium" = "Musculoskeletal", "Connective tissue" = "Musculoskeletal",
  "Sympathetic nerve" = "Neural", "Kidney" = "Visceral organ",
  "Dorsal root ganglion" = "Neural", "GI tract" = "Visceral organ",
  "Cartilage" = "Musculoskeletal", "Mucosal epithelium" = "Epithelial/other",
  "Choroid plexus" = "Neural", "Adipose tissue" = "Epithelial/other",
  "Jaw and tooth" = "Musculoskeletal", "Cavity" = "Epithelial/other",
  "Bone" = "Musculoskeletal", "Heart" = "Cardiopulmonary",
  "Meninges" = "Neural", "Epidermis" = "Epithelial/other"
)

tissue_group_colours <- c(
  "Neural" = "#8EC5D9",
  "Musculoskeletal" = "#F2A65A",
  "Visceral organ" = "#A9D18E",
  "Cardiopulmonary" = "#C7B5E5",
  "Epithelial/other" = "#D9D9D9"
)

trait_group <- c(
  "ADHD" = "Neuropsychiatric", "SCZ" = "Neuropsychiatric",
  "BIP" = "Neuropsychiatric", "MDD" = "Neuropsychiatric",
  "ANX" = "Neuropsychiatric", "IQ" = "Neuropsychiatric",
  "LC" = "Other", "TP" = "Metabolic", "DBP" = "Cardiovascular",
  "CAD" = "Cardiovascular"
)

trait_group_colours <- c(
  "Neuropsychiatric" = "#2878B5",
  "Other" = "#E76F00",
  "Metabolic" = "#2CA02C",
  "Cardiovascular" = "#8E1599"
)

read_spot_annotation <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  spot <- f[["obs/cell_name"]][]
  categories <- f[["obs/annotation/categories"]][]
  codes <- as.integer(f[["obs/annotation/codes"]][])
  tissue <- rep(NA_character_, length(codes))
  valid <- codes >= 0L & codes < length(categories)
  tissue[valid] <- categories[codes[valid] + 1L]
  data.table(spot = spot, tissue = tissue)
}

spot_annotation <- read_spot_annotation(input_h5ad)
stopifnot(nrow(spot_annotation) == 121767L)
stopifnot(setequal(unique(na.omit(spot_annotation$tissue)), tissues))

read_method_result <- function(method, trait) {
  if (method == "SDESE") {
    path <- file.path("sdese", paste0(trait, "_E16.5_gene_marker_score.feather.enrichment.tsv"))
    x <- fread(path, sep = "\t", select = c("Condition", "Adjusted(p)"))
    setnames(x, c("Condition", "Adjusted(p)"), c("spot", "adjusted_p"))
  } else if (method == "SLDSC") {
    path <- file.path("sldsc", paste0("E16.5_", trait, ".csv.gz"))
    x <- fread(path, select = c("spot", "p"))
    # Benjamini-Hochberg FDR correction (equivalent to statsmodels fdr_bh).
    x[, adjusted_p := p.adjust(p, method = "BH")]
    x[, p := NULL]
  } else {
    stop("Unknown method: ", method)
  }
  x[, significant := !is.na(adjusted_p) & adjusted_p < 0.05]
  x
}

compute_enrichment <- function(method, trait) {
  result <- read_method_result(method, trait)

  x <- merge(spot_annotation, result[, .(spot, significant)], by = "spot", all.x = TRUE)
  x[is.na(significant), significant := FALSE]

  total_sig <- sum(x$significant)
  total_n <- nrow(x)
  out <- rbindlist(lapply(tissues, function(tissue_name) {
    in_tissue <- x$tissue == tissue_name
    a <- sum(in_tissue & x$significant, na.rm = TRUE)
    b <- sum(in_tissue & !x$significant, na.rm = TRUE)
    c_count <- total_sig - a
    d <- (total_n - total_sig) - b
    tab <- matrix(c(a, b, c_count, d), nrow = 2, byrow = TRUE)
    fisher_p <- fisher.test(tab, alternative = "greater")$p.value
    odds_ratio <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c_count + 0.5))
    data.table(
      method = method, trait = trait, tissue = tissue_name,
      significant_in_tissue = a, total_in_tissue = a + b,
      significant_outside = c_count, total_outside = c_count + d,
      odds_ratio = odds_ratio, fisher_p = fisher_p,
      total_significant_spots = total_sig,
      selection_rule = ifelse(
        method == "SDESE",
        "Adjusted(p) < 0.05",
        "BH-FDR adjusted p < 0.05"
      )
    )
  }))
  out[, tissue_fdr := p.adjust(fisher_p, method = "BH")]
  # Like the reference, non-enriched/non-significant cells return to the blue baseline.
  out[, plotted_odds_ratio := fifelse(tissue_fdr < 0.05 & odds_ratio > 1, odds_ratio, 0)]
  out[, plotted_odds_ratio_capped := pmin(plotted_odds_ratio, 10)]
  out
}

source_data_path <- file.path(output_dir, "E16.5_tissue_enrichment_source_data.csv")
reuse_source <- identical(Sys.getenv("REUSE_HEATMAP_SOURCE"), "1") && file.exists(source_data_path)
if (reuse_source) {
  message("Reusing existing source-data table: ", source_data_path)
  all_stats <- fread(source_data_path)
} else {
  all_stats <- rbindlist(lapply(c("SDESE", "SLDSC"), function(method) {
    rbindlist(lapply(traits, function(trait) compute_enrichment(method, trait)))
  }))
  fwrite(all_stats, source_data_path)
}

stopifnot(nrow(all_stats) == length(traits) * length(tissues) * 2L)
stopifnot(setequal(unique(all_stats$method), c("SDESE", "SLDSC")))
stopifnot(all(is.finite(all_stats$odds_ratio)))
stopifnot(all(all_stats$tissue_fdr >= 0 & all_stats$tissue_fdr <= 1))

make_heatmap <- function(method_name) {
  d <- copy(all_stats[method == method_name])
  d[, x := match(tissue, tissues)]
  d[, y := length(traits) - match(trait, traits) + 1]

  ann <- data.table(
    tissue = tissues,
    x = seq_along(tissues),
    group = unname(tissue_group[tissues])
  )
  ann[, colour := unname(tissue_group_colours[group])]

  row_ann <- data.table(
    trait = traits,
    y = length(traits) - seq_along(traits) + 1,
    group = unname(trait_group[traits])
  )
  row_ann[, colour := unname(trait_group_colours[group])]

  legend_groups <- names(tissue_group_colours)
  legend_x <- c(1.3, 6.0, 11.2, 16.9, 22.0)

  ggplot(d, aes(x = x, y = y)) +
    geom_tile(aes(fill = plotted_odds_ratio_capped), width = 0.98, height = 0.98,
              colour = "#FFFFFF", linewidth = 0.23) +
    geom_tile(data = ann, aes(x = x, y = 11.05), inherit.aes = FALSE,
              width = 0.98, height = 0.22, fill = ann$colour, colour = NA) +
    geom_tile(data = row_ann, aes(x = 0.28, y = y), inherit.aes = FALSE,
              width = 0.12, height = 0.80, fill = row_ann$colour, colour = NA) +
    geom_text(data = row_ann, aes(x = 0.13, y = y, label = trait), inherit.aes = FALSE,
              hjust = 1, family = "Arial", fontface = "bold", size = 3.15,
              colour = row_ann$colour) +
    annotate("point", x = legend_x, y = 11.72, shape = 15, size = 2.3,
             colour = unname(tissue_group_colours[legend_groups])) +
    annotate("text", x = legend_x + 0.32, y = 11.72, label = legend_groups,
             hjust = 0, family = "Arial", size = 2.35, colour = "#404040") +
    scale_fill_gradientn(
      colours = c("#0B3C6F", "#2C7FB8", "#A6CEE3", "#F7F7F7", "#F4A582", "#B2182B", "#67001F"),
      values = rescale(c(0, 1.5, 3, 4.5, 6, 8, 10)),
      limits = c(0, 10), oob = squish,
      breaks = c(0, 2, 4, 6, 8, 10),
      labels = c("0", "2", "4", "6", "8", "10+"),
      name = "Odds ratio"
    ) +
    scale_x_continuous(
      breaks = seq_along(tissues), labels = tissues,
      limits = c(-1.15, length(tissues) + 0.55), expand = c(0, 0),
      position = "bottom"
    ) +
    scale_y_continuous(limits = c(0.48, 12.03), breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(
      ifelse(method_name == "SDESE", "SDESE Adjusted(p) < 0.05", "S-LDSC BH-FDR < 0.05"),
      " - E16.5 tissue enrichment"
    )) +
    theme_minimal(base_family = "Arial", base_size = 7) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = 6.4, colour = "#222222"),
      axis.ticks.x = element_line(linewidth = 0.25, colour = "#555555"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 6.5),
      legend.key.height = grid::unit(16, "mm"),
      legend.key.width = grid::unit(3.5, "mm"),
      plot.margin = margin(t = 3, r = 4, b = 3, l = 8, unit = "mm")
    ) +
    guides(fill = guide_colourbar(frame.colour = "#555555", frame.linewidth = 0.25,
                                  ticks.colour = "#555555", title.position = "top")) +
    labs(title = NULL)
}

save_heatmap <- function(plot, stem, width_mm = 183, height_mm = 112, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  ragg::agg_png(paste0(stem, ".png"), width = width_in, height = height_in,
                units = "in", res = 300, background = "white")
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
                       family = "Arial", onefile = TRUE)
  print(plot)
  dev.off()
}

for (method in c("SDESE", "SLDSC")) {
  p <- make_heatmap(method)
  save_heatmap(p, file.path(output_dir, paste0("E16.5_", method, "_tissue_enrichment_heatmap")))
}

summary_dt <- all_stats[, .(
  significant_spots = unique(total_significant_spots),
  enriched_tissues_fdr_0.05 = sum(tissue_fdr < 0.05 & odds_ratio > 1),
  maximum_odds_ratio = max(odds_ratio, na.rm = TRUE)
), by = .(method, trait)]
fwrite(summary_dt, file.path(output_dir, "E16.5_heatmap_QA_summary.csv"))

qa_notes <- c(
  "E16.5 tissue-enrichment heatmap QA",
  "Archetype: quantitative grid",
  "Core conclusion: significant spatial localizations show tissue-specific enrichment patterns that differ between SDESE and SLDSC.",
  "Universe: 121,767 annotated spots across 25 tissues.",
  "SDESE significance: Adjusted(p) < 0.05, as in the supplied notebook.",
  "SLDSC significance: Benjamini-Hochberg FDR-adjusted p < 0.05 (fdr_bh).",
  "Cell statistic: one-sided Fisher tissue enrichment; Haldane-Anscombe corrected odds ratio.",
  "Display rule: tissue FDR < 0.05 and odds ratio > 1; otherwise plotted as 0.",
  "Color scale: common 0-10 odds-ratio scale for both methods; values >10 are capped only in the figure.",
  "Uncapped values and exact P/FDR values are retained in E16.5_tissue_enrichment_source_data.csv.",
  "Exports: 183 x 112 mm; PNG preview at 300 dpi,  plus vector PDF.",
  "Visual QA: inspect PNG previews at native size; confirm labels, annotation bars, shared scale, and margins."
)
writeLines(qa_notes, file.path(output_dir, "E16.5_heatmap_QA_notes.txt"), useBytes = TRUE)

cat("Completed heatmaps.\n")
print(summary_dt)
