suppressPackageStartupMessages({
  library(ggplot2)
  library(ggbeeswarm)
  library(dplyr)
  library(readr)
})

base_dir <- Sys.getenv(
  "PRECISION_BASE_DIR",
  unset = "C:\\Users\\Administrator\\Desktop\\SDESE_summary\\simulation"
)
tsv_path <- file.path(base_dir, "metrics_precision.tsv")
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(tsv_path)) {
  stop("Input file does not exist: ", tsv_path)
}

raw_df <- read.delim(
  tsv_path,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

records <- vector("list", nrow(raw_df) * ncol(raw_df))
k <- 0L

for (i in seq_len(nrow(raw_df))) {
  parts <- strsplit(rownames(raw_df)[i], "\\|", fixed = FALSE)[[1]]

  if (length(parts) != 2L) {
    stop("Unexpected row name: ", rownames(raw_df)[i])
  }

  scenario_text <- trimws(parts[1])
  case_text <- trimws(parts[2])

  for (j in seq_len(ncol(raw_df))) {
    cell <- raw_df[i, j]

    if (!is.na(cell) && cell != "NA" && nzchar(cell)) {
      k <- k + 1L
      records[[k]] <- data.frame(
        scenario = scenario_text,
        case = case_text,
        slice = colnames(raw_df)[j],
        Precision = suppressWarnings(as.numeric(cell)),
        stringsAsFactors = FALSE
      )
    }
  }
}

records <- records[seq_len(k)]
long_df <- bind_rows(records) %>%
  filter(is.finite(Precision))

message("Valid precision values: ", nrow(long_df))
message(
  "Precision range: ",
  format(min(long_df$Precision), digits = 4),
  " to ",
  format(max(long_df$Precision), digits = 4)
)

scenario_labels <- c(
  sim_3030 = "1 hotspot",
  sim_1010_9090 = "2 hotspots",
  `sim_r40-59_c40-59` = "1 region",
  `sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99` = "3 regions"
)

scenario_order <- names(scenario_labels)
ratio_order <- c(25L, 50L, 75L)

long_df <- long_df %>%
  mutate(
    scenario_id = sub("\\s+[0-9]+$", "", scenario),
    driver_ratio = as.integer(sub("^.*\\s+([0-9]+)$", "\\1", scenario)),
    scenario_id = factor(scenario_id, levels = scenario_order),
    driver_ratio = factor(driver_ratio, levels = ratio_order),
    case = factor(
      case,
      levels = c("p-value_Condi", "spatial_Condi")
    )
  )

if (any(is.na(long_df$scenario_id))) {
  stop("Some scenarios are not defined in scenario_labels.")
}

if (any(is.na(long_df$driver_ratio))) {
  stop("Some driver ratios are not 25, 50 or 75.")
}

if (any(is.na(long_df$case))) {
  stop("Some case values are not p-value_Condi or spatial_Condi.")
}

case_colours <- c(
  "p-value_Condi" = "#1F5A94",
  "spatial_Condi" = "#D46A4C"
)

case_labels <- c(
  "p-value_Condi" = "P-value condition",
  "spatial_Condi" = "Spatial condition"
)

median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(data.frame(y = NA_real_, ymin = NA_real_, ymax = NA_real_))
  }

  data.frame(
    y = median(x),
    ymin = unname(quantile(x, 0.25)),
    ymax = unname(quantile(x, 0.75))
  )
}

set.seed(20260831)

precision_plot <- ggplot(
  long_df,
  aes(
    x = driver_ratio,
    y = Precision,
    colour = case,
    fill = case,
    group = case
  )
) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 0.52,
    width = 0.10,
    size = 0.58,
    alpha = 0.24,
    stroke = 0,
    na.rm = TRUE
  ) +
  stat_summary(
    fun.data = median_iqr,
    geom = "linerange",
    position = position_dodge(width = 0.52),
    linewidth = 0.80,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.52),
    shape = 21,
    size = 2.0,
    stroke = 0.35,
    colour = "white",
    na.rm = TRUE
  ) +
  facet_grid(
    cols = vars(scenario_id),
    labeller = labeller(scenario_id = as_labeller(scenario_labels))
  ) +
  scale_colour_manual(
    values = case_colours,
    labels = case_labels,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = case_colours,
    labels = case_labels,
    drop = FALSE
  ) +
  scale_x_discrete(labels = c("25", "50", "75")) +
  scale_y_continuous(
    breaks = c(0.25, 0.50, 0.75, 1.00),
    labels = c("0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  coord_cartesian(ylim = c(0.25, 1.02)) +
  labs(
    x = "Driver ratio (%)",
    y = "Precision",
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(colour = "#333333", linewidth = 0.35),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.2, "mm"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6.3, colour = "#333333"),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    strip.text = element_text(size = 6.8, face = "bold"),
    panel.spacing.x = grid::unit(0.9, "mm"),
    panel.grid.major.y = element_line(colour = "#E7E7E7", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(6, "mm"),
    plot.margin = margin(3, 4, 3, 4)
  )

width_in <- 112 / 25.4
height_in <- 62 / 25.4

ggsave(
  file.path(fig_dir, "genes_precision_raw_median_iqr.png"),
  precision_plot,
  width = width_in,
  height = height_in,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "genes_precision_raw_median_iqr.pdf"),
  precision_plot,
  width = width_in,
  height = height_in,
  bg = "white",
  device = grDevices::cairo_pdf,
  family = "Arial"
)

write_csv(
  long_df %>%
    mutate(
      scenario_id = as.character(scenario_id),
      driver_ratio = as.integer(as.character(driver_ratio)),
      case = as.character(case)
    ) %>%
    arrange(scenario_id, driver_ratio, case, slice),
  file.path(fig_dir, "genes_precision_source_data.csv")
)

message("Saved precision figure to: ", fig_dir)
