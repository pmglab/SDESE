suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggbeeswarm)
})

# Paths and export settings ---------------------------------------------------
base_dir <- Sys.getenv(
  "SDESE_BASE_DIR",
  unset = "C:\\Users\\Administrator\\Desktop\\SDESE_summary\\simulation"
)

input_file <- file.path(base_dir, "sensitivity_nsr_summary.tsv")
output_dir <- file.path(base_dir, "figures", "sensitivity_nsr_summary")
details_dir <- file.path(output_dir, "details")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(details_dir, recursive = TRUE, showWarnings = FALSE)

total_spots <- 10000
png_dpi <- 600
save_pdf <- TRUE

method_colours <- c(SDESE = "#1F5A94", SLDSC = "#D46A4C")
method_shapes <- c(SDESE = 16, SLDSC = 17)
method_lines <- c(SDESE = "solid", SLDSC = "22")

red_colours <- c(
  "#FFF7F3", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
  "#EF3B2C", "#CB181D", "#A50F15", "#67000D"
)

nsr_colours <- c(
  "#053061", "#2166AC", "#67A9CF",
  "#D9D9D9",
  "#EF8A62", "#B2182B", "#67001F"
)

nsr_colour_values <- scales::rescale(
  c(0, 0.20, 0.40, 0.50, 0.65, 0.82, 1.00)
)

scenario_levels <- c(
  "sim_3030",
  "sim_1010_9090",
  "sim_r40-59_c40-59",
  "sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99"
)

scenario_labels <- c(
  sim_3030 = "1 hotspot",
  sim_1010_9090 = "2 hotspots",
  `sim_r40-59_c40-59` = "1 region",
  `sim_r20-39_c20-39_r50-79_c50-79_r90-99_c90-99` = "3 regions"
)

group_levels <- unlist(lapply(
  scenario_levels,
  function(x) paste0(scenario_labels[[x]], "\n", c(25, 50, 75), "%")
))

data <- read_tsv(
  input_file,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  filter(toupper(method) %in% c("SDESE", "SLDSC")) %>%
  mutate(
    method = factor(toupper(method), levels = c("SDESE", "SLDSC")),
    scenario = factor(scenario, levels = scenario_levels),
    driver_ratio = as.integer(driver_ratio),
    mean_sensitivity = as.numeric(mean_sensitivity),
    mean_NSR = as.numeric(mean_NSR),
    significant_spot_simulation_rate =
      as.numeric(significant_spot_simulation_rate)
  )

required <- c(
  "method", "scenario", "driver_ratio", "mean_sensitivity", "mean_NSR",
  "significant_spot_simulation_rate", "simulation_significant_spot_counts",
  "simulation_NSR_values"
)

missing_columns <- setdiff(required, names(data))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

if (any(is.na(data$scenario))) {
  stop("The input contains scenarios not defined in scenario_levels.")
}

add_group <- function(x) {
  x %>%
    mutate(
      group = paste0(
        scenario_labels[as.character(scenario)],
        "\n",
        driver_ratio,
        "%"
      ),
      group = factor(group, levels = group_levels)
    )
}

expand_values <- function(data, column, value_name) {
  rows <- lapply(seq_len(nrow(data)), function(i) {
    cell <- data[[column]][i]

    if (is.na(cell) || !nzchar(cell)) {
      stop(
        data$method[i], " / ", data$scenario[i], " / ",
        data$driver_ratio[i], "%: ", column, " is empty."
      )
    }

    values <- strsplit(cell, ";", fixed = TRUE)[[1]]
    if (length(values) != 100L) {
      stop(
        data$method[i], " / ", data$scenario[i], " / ",
        data$driver_ratio[i], "%: ", column,
        " must contain exactly 100 values."
      )
    }

    values[values %in% c("", "NA", "NaN", "nan")] <- NA_character_

    tibble(
      method = data$method[i],
      scenario = data$scenario[i],
      driver_ratio = data$driver_ratio[i],
      simulation = 0:99,
      value = suppressWarnings(as.numeric(values))
    )
  })

  result <- bind_rows(rows)
  names(result)[names(result) == "value"] <- value_name
  add_group(result)
}

count_long <- expand_values(
  data,
  "simulation_significant_spot_counts",
  "spot_count"
)

nsr_long <- expand_values(
  data,
  "simulation_NSR_values",
  "NSR"
)

summary_data <- add_group(data)

save_plot <- function(
    plot,
    name,
    width = 8,
    height = 7,
    folder = output_dir) {
  ggsave(
    file.path(folder, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = png_dpi,
    bg = "white"
  )

  if (save_pdf) {
    ggsave(
      file.path(folder, paste0(name, ".pdf")),
      plot,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      device = grDevices::cairo_pdf,
      family = "Arial"
    )
  }
}

# Distribution plots
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

distribution_theme <- theme_classic(base_size = 7, base_family = "Arial") +
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
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5),
    legend.key.width = grid::unit(6, "mm"),
    plot.margin = margin(3, 4, 3, 4)
  )

set.seed(20260831)

count_distribution_plot <- ggplot(
  count_long,
  aes(
    x = factor(driver_ratio, levels = c(25, 50, 75)),
    y = spot_count,
    colour = method,
    fill = method,
    group = method
  )
) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 0.52,
    width = 0.10,
    size = 0.62,
    alpha = 0.24,
    stroke = 0,
    na.rm = TRUE
  ) +
  stat_summary(
    fun.data = median_iqr,
    geom = "linerange",
    position = position_dodge(width = 0.52),
    linewidth = 0.85,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.52),
    shape = 21,
    size = 2.15,
    stroke = 0.35,
    colour = "white",
    na.rm = TRUE
  ) +
  facet_grid(
    cols = vars(scenario),
    labeller = labeller(scenario = as_labeller(scenario_labels))
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_fill_manual(values = method_colours, drop = FALSE) +
  scale_x_discrete(labels = c("25", "50", "75")) +
  scale_y_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 10, 100, 1000, 2500),
    labels = c("0", "10", "100", "1,000", "2,500"),
    expand = expansion(mult = c(0.03, 0.07))
  ) +
  labs(
    x = "Driver ratio (%)",
    y = "Significant spots per slice\n(log1p scale)"
  ) +
  distribution_theme

nsr_distribution_plot <- ggplot(
  nsr_long,
  aes(
    x = factor(driver_ratio, levels = c(25, 50, 75)),
    y = NSR,
    colour = method,
    fill = method,
    group = method
  )
) +
  ggbeeswarm::geom_quasirandom(
    dodge.width = 0.52,
    width = 0.10,
    size = 0.62,
    alpha = 0.24,
    stroke = 0,
    na.rm = TRUE
  ) +
  stat_summary(
    fun.data = median_iqr,
    geom = "linerange",
    position = position_dodge(width = 0.52),
    linewidth = 0.85,
    na.rm = TRUE
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.52),
    shape = 21,
    size = 2.15,
    stroke = 0.35,
    colour = "white",
    na.rm = TRUE
  ) +
  facet_grid(
    cols = vars(scenario),
    labeller = labeller(scenario = as_labeller(scenario_labels))
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_fill_manual(values = method_colours, drop = FALSE) +
  scale_x_discrete(labels = c("25", "50", "75")) +
  scale_y_continuous(
    limits = c(-0.04, 1.05),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Driver ratio (%)", y = "NSR") +
  distribution_theme

write_csv(
  bind_rows(
    count_long %>%
      transmute(
        method,
        scenario,
        driver_ratio,
        simulation,
        metric = "Significant spot count",
        value = spot_count
      ),
    nsr_long %>%
      transmute(
        method,
        scenario,
        driver_ratio,
        simulation,
        metric = "NSR",
        value = NSR
      )
  ),
  file.path(output_dir, "simulation_distribution_source_data.csv")
)

decimal2_label <- function(x) formatC(x, format = "f", digits = 2)

rate_data <- summary_data %>%
  select(
    method,
    scenario,
    driver_ratio,
    significant_spot_simulation_rate
  ) %>%
  rename(value = significant_spot_simulation_rate) %>%
  mutate(
    value_label = decimal2_label(value),
    label_vjust = case_when(
      method == "SDESE" ~ 1.65,
      method == "SLDSC" & value >= 0.80 ~ 3.15,
      TRUE ~ -0.85
    )
  )

significant_rate_plot <- ggplot(
  rate_data,
  aes(
    x = driver_ratio,
    y = value,
    colour = method,
    shape = method,
    linetype = method,
    group = method
  )
) +
  geom_hline(
    yintercept = c(0.25, 0.50, 0.75),
    colour = "#E7E7E7",
    linewidth = 0.25
  ) +
  geom_line(linewidth = 0.65, na.rm = TRUE) +
  geom_point(size = 2.15, stroke = 0.25, na.rm = TRUE) +
  geom_text(
    aes(label = value_label, vjust = label_vjust),
    size = 1.80,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  facet_grid(
    cols = vars(scenario),
    labeller = labeller(scenario = as_labeller(scenario_labels))
  ) +
  scale_x_continuous(
    breaks = c(25, 50, 75),
    labels = c("25", "50", "75"),
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  scale_y_continuous(
    limits = c(-0.06, 1.08),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_colour_manual(values = method_colours, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  scale_linetype_manual(values = method_lines, drop = FALSE) +
  labs(
    x = "Driver ratio (%)",
    y = "Proportion of slices\nwith significant spots",
    colour = NULL,
    shape = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(colour = "#333333", linewidth = 0.35),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.3, "mm"),
    axis.text = element_text(colour = "#333333", size = 6.5),
    axis.text.y = element_text(margin = margin(r = 2.5)),
    axis.title.x = element_text(size = 7, margin = margin(t = 5)),
    strip.background.x = element_rect(fill = "#F2F2F2", colour = NA),
    strip.text.x = element_text(
      face = "bold",
      size = 7,
      margin = margin(3, 2, 3, 2)
    ),
    panel.spacing.x = grid::unit(1.4, "mm"),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 6.8),
    legend.key.width = grid::unit(7, "mm"),
    legend.spacing.x = grid::unit(2, "mm"),
    plot.margin = margin(4, 5, 4, 4)
  )

sensitivity_data <- summary_data %>%
  select(method, scenario, driver_ratio, mean_sensitivity) %>%
  rename(value = mean_sensitivity) %>%
  mutate(
    value_label = decimal2_label(value),
    label_vjust = case_when(
      method == "SDESE" ~ 1.65,
      method == "SLDSC" & value >= 0.80 ~ 3.15,
      TRUE ~ -0.85
    )
  )

sensitivity_line_plot <- significant_rate_plot %+% sensitivity_data +
  labs(y = "sensitivity")

save_plot(
  significant_rate_plot,
  "significant_spot_simulation_rate",
  width = 112 / 25.4,
  height = 66 / 25.4
)

save_plot(
  count_distribution_plot,
  "simulation_significant_spot_counts",
  width = 112 / 25.4,
  height = 66 / 25.4
)

save_plot(
  sensitivity_line_plot,
  "simulation_sensitivity_values",
  width = 112 / 25.4,
  height = 66 / 25.4
)

save_plot(
  nsr_distribution_plot,
  "simulation_NSR_values",
  width = 112 / 25.4,
  height = 66 / 25.4
)

write_csv(
  rate_data %>%
    mutate(scenario = as.character(scenario)) %>%
    arrange(scenario, driver_ratio, method),
  file.path(output_dir, "significant_rate_source_data.csv")
)

write_csv(
  sensitivity_data %>%
    mutate(scenario = as.character(scenario)) %>%
    arrange(scenario, driver_ratio, method),
  file.path(output_dir, "sensitivity_source_data.csv")
)

heatmap_group_levels <- gsub("\n", " · ", group_levels, fixed = TRUE)

heatmap_rows <- unlist(lapply(
  heatmap_group_levels,
  function(x) paste(x, c("SDESE", "SLDSC"), sep = " | ")
))

prepare_heatmap <- function(data) {
  data %>%
    mutate(
      heatmap_group = gsub("\n", " · ", as.character(group), fixed = TRUE),
      row = factor(
        paste(heatmap_group, method, sep = " | "),
        levels = rev(heatmap_rows)
      )
    )
}

percentage_heatmap <- function(data, scenarios) {
  plot_data <- data %>%
    filter(scenario %in% scenarios) %>%
    mutate(percentage = spot_count / total_spots * 100) %>%
    prepare_heatmap()

  limits <- range(plot_data$percentage, na.rm = TRUE)
  if (diff(limits) == 0) {
    limits[2] <- limits[1] + 1e-6
  }

  n_scenarios <- dplyr::n_distinct(plot_data$scenario)
  n_heatmap_rows <- n_scenarios * 6L

  pair_separators <- seq(
    from = 2.5,
    to = n_heatmap_rows - 0.5,
    by = 2
  )

  scenario_separators <- seq(
    from = 6.5,
    to = n_heatmap_rows - 0.5,
    by = 6
  )

  ratio_separators <- setdiff(pair_separators, scenario_separators)

  ggplot(plot_data, aes(simulation, row, fill = percentage)) +
    geom_tile() +
    geom_hline(
      yintercept = ratio_separators,
      colour = "white",
      linewidth = 4.0
    ) +
    geom_hline(
      yintercept = scenario_separators,
      colour = "white",
      linewidth = 4.0
    ) +
    scale_x_continuous(
      limits = c(-0.5, 99.5),
      breaks = c(0, 20, 40, 60, 80, 99),
      expand = c(0, 0)
    ) +
    scale_fill_gradientn(
      colours = red_colours,
      limits = limits,
      labels = function(x) paste0(format(x, digits = 3), "%"),
      name = "Significant\nspots"
    ) +
    labs(x = "Simulated slice ID", y = NULL) +
    theme_minimal(base_size = 14, base_family = "Helvetica") +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
}

percentage_plot <- percentage_heatmap(
  count_long,
  scenario_levels
)

# Complete four-scenario heatmap, matching the NSR heatmap layout.
save_plot(
  percentage_plot,
  "simulation_significant_spot_percentage",
  width = 16,
  height = 8,
  folder = details_dir
)

nsr_heatmap_data <- prepare_heatmap(nsr_long)

nsr_pair_separators <- seq(2.5, 23.5, by = 2)
nsr_scenario_separators <- c(6.5, 12.5, 18.5)
nsr_ratio_separators <- setdiff(
  nsr_pair_separators,
  nsr_scenario_separators
)

nsr_heatmap <- ggplot(
  nsr_heatmap_data,
  aes(simulation, row, fill = NSR)
) +
  geom_tile() +
  geom_hline(
    yintercept = nsr_ratio_separators,
    colour = "white",
    linewidth = 4.0
  ) +
  geom_hline(
    yintercept = nsr_scenario_separators,
    colour = "white",
    linewidth = 4.0
  ) +
  scale_x_continuous(
    limits = c(-0.5, 99.5),
    breaks = c(0, 20, 40, 60, 80, 99),
    expand = c(0, 0)
  ) +
  scale_fill_gradientn(
    colours = nsr_colours,
    values = nsr_colour_values,
    limits = c(0, 1),
    na.value = "#A6A6A6",
    breaks = seq(0, 1, 0.25),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    oob = scales::squish,
    name = "NSR"
  ) +
  labs(x = "Simulated slice ID", y = NULL) +
  theme_minimal(base_size = 14, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

save_plot(
  nsr_heatmap,
  "simulation_NSR_values",
  width = 16,
  height = 8,
  folder = details_dir
)

message("All figures saved to: ", output_dir)
