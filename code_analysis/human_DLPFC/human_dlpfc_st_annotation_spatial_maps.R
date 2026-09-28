suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(hdf5r)
  library(jsonlite)
  library(png)
})

project_root <- normalizePath(getwd(), winslash = "/")
h5ad_dir <- file.path(project_root, "resources", "Brain_33558695_h5ad")
out_dir <- file.path(project_root, "plot", "human_dlpfc_st_annotation_spatial_maps")
ffmpeg <- "C:/Program Files/ffmpeg/bin/ffmpeg.exe"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- paste0("V", 1:12)
files <- file.path(h5ad_dir, sprintf("HS_Brain_33558695_%s.seuratobj.h5ad", samples))
image_files <- file.path(h5ad_dir, sprintf("HS_Brain_33558695_%s.spatial_image.jpg", samples))
scale_files <- file.path(h5ad_dir, sprintf("HS_Brain_33558695_%s.imagescale.json", samples))
hires_scales <- map_dbl(scale_files, ~ fromJSON(.x)$hires[[1]])

read_obs <- function(h5, column) {
  x <- h5[[paste0("obs/", column)]]
  if (inherits(x, "H5Group")) {
    codes <- as.integer(x[["codes"]]$read())
    categories <- as.character(x[["categories"]]$read())
    return(ifelse(codes < 0, NA_character_, categories[codes + 1L]))
  }
  x$read()
}

spot_data <- map_dfr(seq_along(samples), function(i) {
  path <- files[i]
  sample_id <- samples[i]
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  tibble(
    sample = sample_id,
    spot = as.character(read_obs(h5, "_index")),
    x = as.numeric(read_obs(h5, "coord_x")),
    y = as.numeric(read_obs(h5, "coord_y")),
    hires_scale = hires_scales[i],
    image_row = x * hires_scale,
    image_col = y * hires_scale,
    spatial_region = paste0("SR", read_obs(h5, "spatial_region")),
    cell_type = as.character(read_obs(h5, "cell_type"))
  )
}) %>%
  mutate(sample = factor(sample, levels = samples))

write_tsv(spot_data, file.path(out_dir, "human_dlpfc_st_annotation_spatial_maps.tsv"), na = "NA")
walk(samples, function(sample_id) {
  write_tsv(
    filter(spot_data, sample == sample_id),
    file.path(out_dir, paste0(sample_id, "_st_annotations.tsv")), na = "NA"
  )
})

region_levels <- paste0("SR", 1:7)
region_colors <- setNames(c("#3B4CC0", "#5D7CE6", "#86B6EB", "#B8DDD5", "#F1D57A", "#E88B51", "#B40426"), region_levels)

read_tissue_image <- function(path) {
  temp <- tempfile(fileext = ".png")
  system2(ffmpeg, c("-y", "-loglevel", "quiet", "-i", shQuote(path), shQuote(temp)))
  image <- readPNG(temp)
  unlink(temp)
  image
}

make_map <- function(data, image) {
  width <- dim(image)[2]
  height <- dim(image)[1]
  data <- mutate(data, plot_x = image_col, plot_y = height - image_row)
  ggplot(data, aes(plot_x, plot_y, color = spatial_region)) +
    annotation_raster(image, xmin = 0, xmax = width, ymin = 0, ymax = height, interpolate = TRUE) +
    geom_point(size = 1.7, alpha = 0.78, stroke = 0) +
    scale_color_manual(values = region_colors, na.value = "#D9D9D9", drop = FALSE) +
    coord_fixed(xlim = c(0, width), ylim = c(0, height), expand = FALSE) +
    labs(title = as.character(unique(data$sample)), color = "Spatial region") +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key.height = grid::unit(4.2, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      legend.position = "right",
      plot.margin = margin(8, 8, 8, 8)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3.2), ncol = 1))
}

walk(samples, function(sample_id) {
  data <- filter(spot_data, sample == sample_id)
  image <- read_tissue_image(image_files[match(sample_id, samples)])
  p <- make_map(data, image)
  stem <- file.path(out_dir, paste0(sample_id, "_spatial_region_spatial_map"))
  ggsave(paste0(stem, ".png"), p, width = 7.2, height = 6.2, dpi = 600, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = 7.2, height = 6.2,
         device = grDevices::cairo_pdf, bg = "white")
})
