suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(hdf5r)
  library(ragg)
  library(scales)
  library(scattermore)
})

options(stringsAsFactors = FALSE)

input_h5ad <- "E16.5_E1S1.MOSTA.h5ad"
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "figures")
all_traits <- c("ADHD", "SCZ", "BIP", "MDD", "ANX", "IQ", "LC", "TP", "DBP", "CAD")
all_methods <- c("SDESE", "SLDSC")
traits <- if (nzchar(Sys.getenv("TRAITS"))) strsplit(Sys.getenv("TRAITS"), ",", fixed = TRUE)[[1]] else all_traits
methods <- if (nzchar(Sys.getenv("METHODS"))) strsplit(Sys.getenv("METHODS"), ",", fixed = TRUE)[[1]] else all_methods
stopifnot(all(traits %in% all_traits), all(methods %in% all_methods))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Common visual contract across all traits and both methods.
background_colour <- "#F0F8FF"
all_spot_colour <- "#6495ED"
signal_colours <- c("#FFD700", "#FF8C00", "#FF4500", "#B22222", "#6E001A")
lower_colour_limit <- -log10(0.05)
upper_colour_limit <- 5
plot_width_mm <- 90
plot_height_mm <- 145

read_spatial_data <- function(path) {
  f <- H5File$new(path, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  spot <- f[["obs/cell_name"]][]
  spatial <- f[["obsm/spatial"]][, ]
  if (nrow(spatial) != length(spot) && ncol(spatial) == length(spot)) {
    spatial <- t(spatial)
  }
  stopifnot(nrow(spatial) == length(spot), ncol(spatial) >= 2L)

  data.table(
    spot = as.character(spot),
    x = as.numeric(spatial[, 1]),
    y = as.numeric(spatial[, 2])
  )
}

read_adjusted_p <- function(method, trait) {
  if (method == "SDESE") {
    path <- file.path(
      "sdese",
      paste0(trait, "_E16.5_gene_marker_score.feather.enrichment.tsv")
    )
    # These files end with a plain-text method note; fread correctly discards it.
    result <- suppressWarnings(
      fread(path, sep = "\t", select = c("Condition", "Adjusted(p)"))
    )
    setnames(result, c("Condition", "Adjusted(p)"), c("spot", "adjusted_p"))
  } else if (method == "SLDSC") {
    path <- file.path("sldsc", paste0("E16.5_", trait, ".csv.gz"))
    result <- fread(path, select = c("spot", "p"))
    result[, adjusted_p := p.adjust(p, method = "BH")]
    result[, p := NULL]
  } else {
    stop("Unknown method: ", method)
  }

  result[, spot := as.character(spot)]
  result[, adjusted_p := pmin(pmax(as.numeric(adjusted_p), 0), 1)]
  result
}

make_spatial_plot <- function(spatial_data, method, trait) {
  result <- read_adjusted_p(method, trait)
  plot_data <- merge(spatial_data, result, by = "spot", all.x = TRUE, sort = FALSE)
  plot_data[, significant := is.finite(adjusted_p) & adjusted_p < 0.05]
  plot_data[, neg_log10_adjusted_p := -log10(pmax(adjusted_p, .Machine$double.xmin))]
  plot_data[, colour_value := pmin(neg_log10_adjusted_p, upper_colour_limit)]

  significant_data <- plot_data[significant == TRUE]
  setorder(significant_data, colour_value)
  significant_count <- nrow(significant_data)

  method_label <- method

  p <- ggplot() +
    geom_scattermore(
      data = plot_data,
      aes(x = x, y = y),
      colour = all_spot_colour,
      alpha = 0.42,
      pointsize = 1.35,
      pixels = c(1100, 1800)
    )

  if (significant_count > 0L) {
    p <- p +
      geom_scattermore(
        data = significant_data,
        aes(x = x, y = y, colour = colour_value),
        alpha = 0.90,
        pointsize = 1.55,
        pixels = c(1100, 1800)
      ) +
      scale_colour_gradientn(
        colours = signal_colours,
        limits = c(lower_colour_limit, upper_colour_limit),
        oob = squish,
        breaks = c(lower_colour_limit, 2, 3, 4, 5),
        labels = c("1.30", "2", "3", "4", "5+"),
        name = expression(-log[10](adjusted~p))
      )
  }

  p +
    coord_fixed(expand = FALSE) +
    labs(
      title = paste0(trait, " — ", method_label),
      subtitle = paste0(format(significant_count, big.mark = ","), " significant spots")
    ) +
    annotate(
      "label",
      x = -Inf, y = Inf,
      label = paste0(format(significant_count, big.mark = ","), " significant spots"),
      hjust = -0.08, vjust = 1.25,
      size = 2.1,
      family = "Arial",
      colour = "#4D4D4D",
      fill = alpha("white", 0.78),
      linewidth = 0,
      label.padding = grid::unit(1.0, "mm")
    ) +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      panel.background = element_rect(fill = background_colour, colour = NA),
      plot.background = element_rect(fill = background_colour, colour = NA),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      legend.key.height = grid::unit(10, "mm"),
      legend.key.width = grid::unit(2.8, "mm"),
      plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "mm")
    ) +
    guides(colour = guide_colourbar(
      title.position = "top",
      frame.colour = "#555555",
      frame.linewidth = 0.25,
      ticks.colour = "#555555"
    )) +
    labs(title = NULL, subtitle = NULL)
}

save_spatial_plot <- function(plot, stem) {
  width_in <- plot_width_mm / 25.4
  height_in <- plot_height_mm / 25.4

  ragg::agg_png(
    paste0(stem, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = 600,
    background = background_colour
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "Arial",
    onefile = TRUE,
    bg = background_colour
  )
  print(plot)
  dev.off()
}

spatial_data <- read_spatial_data(input_h5ad)
stopifnot(nrow(spatial_data) == 121767L)

summary_rows <- list()
index <- 1L
for (trait in traits) {
  for (method in methods) {
    message("Plotting ", trait, " / ", method)
    adjusted <- read_adjusted_p(method, trait)
    significant_count <- adjusted[is.finite(adjusted_p) & adjusted_p < 0.05, .N]
    summary_rows[[index]] <- data.table(
      trait = trait,
      method = method,
      significant_spots = significant_count
    )
    index <- index + 1L

    p <- make_spatial_plot(spatial_data, method, trait)
    stem <- file.path(output_dir, paste0(trait, "_", method, "_spatial_visualization"))
    save_spatial_plot(p, stem)
  }
}

summary_table <- rbindlist(summary_rows)
stopifnot(nrow(summary_table) == length(traits) * length(methods))
cat("Completed ", nrow(summary_table), " spatial maps (PNG + PDF).\n", sep = "")
print(summary_table)
