packages <- c("MASS", "ggplot2", "dplyr", "CompQuadForm", "Matrix", "nnls", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("/public3/ly/SDESE_summary/aECS")

create_AR1 <- function(size, rho) {
  outer(1:size, 1:size, function(i, j) rho^abs(i - j))
}

calc_liu_p <- function(q, c1, c2, c3) {
  if (c2 <= 1e-8) return(ifelse(q > c1, 0, 1))
  if (abs(c3) < 1e-8) return(pnorm(q, mean = c1, sd = sqrt(2 * c2), lower.tail = FALSE))

  if (c3 > 0) {
    s <- 2 * c3 / c2
    a <- (c2^3) / (2 * c3^2)
    delta <- c1 - a * s
    pval <- pgamma(q - delta, shape = a, scale = s, lower.tail = FALSE)
  } else {
    s <- 2 * (-c3) / c2
    a <- (c2^3) / (2 * (-c3)^2)
    delta <- c1 + a * s
    pval <- pgamma(delta - q, shape = a, scale = s, lower.tail = TRUE)
  }
  return(max(min(pval, 1), 0))
}

get_beta_mat_from_h2_strict <- function(nSnp, h2, nSim, R_true, causal_frac = 0.2) {
  if (h2 <= 0) return(matrix(0, nrow = nSnp, ncol = nSim))
  n_causal <- max(1, floor(nSnp * causal_frac))
  b_mat <- matrix(0, nrow = nSnp, ncol = nSim)
  for (i in 1:nSim) {
    causal_idx <- sample(1:nSnp, n_causal)
    b_raw <- rep(0, nSnp)
    b_raw[causal_idx] <- rnorm(n_causal)
    v_g <- as.numeric(t(b_raw) %*% R_true %*% b_raw)
    if (v_g > 0) {
      b_raw <- b_raw * sqrt(h2 / v_g)
    }
    b_mat[, i] <- b_raw
  }
  return(b_mat)
}

run_sim_raw_p <- function(nSim = 1000, nSnpA = 50, nSnpB = 30, N = 50000,
                          h2_A = 0.01, h2_B = 0.0, rho_true = 0.4,
                          aECS_c = 1, sign_mismatch_prop = 0.0, ld_noise_sd = 0.0) {

  totalSnp <- nSnpA + nSnpB

  # 1. Generate the true LD matrix
  R_true <- create_AR1(totalSnp, rho_true)
  L_true <- t(chol(R_true))
  R_AA_true <- R_true[1:nSnpA, 1:nSnpA, drop=FALSE]
  R_BB_true <- R_true[(nSnpA+1):totalSnp, (nSnpA+1):totalSnp, drop=FALSE]

  # 2. Inject LD noise and phase-flip mismatch
  R_ref_raw <- R_true
  if (sign_mismatch_prop > 0) {
    flip_mask <- matrix(1, totalSnp, totalSnp)
    upper_idx <- which(upper.tri(flip_mask))
    n_flips <- floor(length(upper_idx) * sign_mismatch_prop)
    if (n_flips > 0) {
      flip_chosen <- sample(upper_idx, n_flips)
      flip_mask[flip_chosen] <- -1
      flip_mask[lower.tri(flip_mask)] <- t(flip_mask)[lower.tri(flip_mask)]
    }
    R_ref_raw <- R_ref_raw * flip_mask
  }

  if (ld_noise_sd > 0) {
    noise_upper <- matrix(0, totalSnp, totalSnp)
    noise_upper[upper.tri(noise_upper)] <- rnorm(sum(upper.tri(noise_upper)), mean = 0, sd = ld_noise_sd)
    noise_mat <- noise_upper + t(noise_upper)
    R_ref_raw <- R_ref_raw + noise_mat
  }

  # 3. Core branch: data routing for Z-score and aECS
  if (sign_mismatch_prop > 0 || ld_noise_sd > 0) {
    # Constrain boundaries
    R_ref_raw[R_ref_raw > 1] <- 1
    R_ref_raw[R_ref_raw < -1] <- -1
    diag(R_ref_raw) <- 1

    # [Branch A]: Traditional Z-score is extremely fragile, so the data must be forcibly distorted to ensure positive semidefiniteness (PSD)
    R_ref_psd <- as.matrix(Matrix::nearPD(R_ref_raw, corr = TRUE, ensureSymmetry = TRUE)$mat)
  } else {
    R_ref_raw <- R_ref_raw
    R_ref_psd <- R_ref_raw
  }

  # Control group: traditional Z-score test (based on R_ref_psd distorted by nearPD)
  R_AA_ref <- R_ref_psd[1:nSnpA, 1:nSnpA, drop=FALSE]
  R_BB_ref <- R_ref_psd[(nSnpA+1):totalSnp, (nSnpA+1):totalSnp, drop=FALSE]
  R_AB_ref <- R_ref_psd[1:nSnpA, (nSnpA+1):totalSnp, drop=FALSE]
  R_BA_ref <- t(R_AB_ref)

  ridgeDelta <- 0.05
  R_AA_ridge <- R_AA_ref + diag(ridgeDelta, nSnpA)
  R_AA_inv_ridge <- solve(R_AA_ridge)
  W_ridge <- R_BA_ref %*% R_AA_inv_ridge

  Sigma_cond_ridge <- R_BB_ref - 2 * (W_ridge %*% R_AB_ref) + W_ridge %*% R_AA_ref %*% t(W_ridge)
  Sigma_cond_ridge_sym <- 0.5 * (Sigma_cond_ridge + t(Sigma_cond_ridge))
  eig_res <- eigen(Sigma_cond_ridge_sym, symmetric = TRUE)
  lambdaArr <- eig_res$values[eig_res$values > 1e-5]

  # Generate the true signal (based on True LD)
  Z_raw <- matrix(rnorm(totalSnp * nSim), nrow = totalSnp, ncol = nSim)
  Z_noise <- L_true %*% Z_raw

  beta_A_mat <- get_beta_mat_from_h2_strict(nSnpA, h2_A, nSim, R_AA_true)
  beta_B_mat <- get_beta_mat_from_h2_strict(nSnpB, h2_B, nSim, R_BB_true)
  beta_full_mat <- rbind(beta_A_mat, beta_B_mat)

  E_Z_mat <- sqrt(N) * (R_true %*% beta_full_mat)
  Z_total <- Z_noise + E_Z_mat

  Z_A <- Z_total[1:nSnpA, , drop=FALSE]
  Z_B <- Z_total[(nSnpA+1):totalSnp, , drop=FALSE]

  # Calculate Z-Score P-values
  mu_ridge <- W_ridge %*% Z_A
  Z_cond_ridge <- Z_B - mu_ridge
  Q_cond_ridge <- colSums(Z_cond_ridge^2)
  pCond_Ridge <- sapply(Q_cond_ridge, function(q) {
    res <- CompQuadForm::davies(q, lambdaArr)
    pval <- res$Qq
    if (res$ifault != 0 || pval <= 0 || pval >= 1) pval <- CompQuadForm::liu(q, lambdaArr)
    return(max(min(pval, 1), 0))
  })

  get_W <- function(R_mat, c_val) {
    n <- nrow(R_mat)
    ridge_eps <- if(c_val < 1) 0.02 else 0.01
    R_stable <- R_mat * (1 - ridge_eps) + diag(ridge_eps, n)
    nnls::nnls(abs(R_stable)^c_val, rep(1, n))$x
  }

  # Independently extract the region A block from the Raw matrix for aECS
  R_AA_raw <- R_ref_raw[1:nSnpA, 1:nSnpA, drop=FALSE]

  W_full <- get_W(R_ref_raw, aECS_c)
  W_A <- get_W(R_AA_raw, aECS_c)

  # 1. Construct differential weights containing negative values (strictly no truncation!)
  W_diff <- W_full - c(W_A, rep(0, nSnpB))

  # 2. Pure Z^2-space mapping to extract the conditional statistic
  S_cond_aECS <- as.numeric(t(W_diff) %*% (Z_total^2))

  # 3. [Absolute-value topology third-moment calibration]
  # Force the use of abs(R_ref_raw) to calculate moments, completely eliminating phase sensitivity,
  # c2 (variance) remains perfectly mathematically equivalent, while c3 (skewness) undergoes forced topological smoothing.
  M <- W_diff * abs(R_ref_raw)
  c1 <- sum(diag(M))
  c2 <- sum(M * t(M))
  M2 <- M %*% M
  c3 <- sum(M2 * t(M))

  # 4. Map to Liu's distribution supporting indefinite quadratic forms (3-Moment Shifted Gamma)
  pCond_aECS <- sapply(S_cond_aECS, function(q) calc_liu_p(q, c1, c2, c3))

  # Consolidate and return
  res_raw <- data.frame(
    SimID = 1:nSim,
    nSnpA = nSnpA,
    Rho = rho_true,
    Type = ifelse(h2_B == 0, "Null Hypothesis", "Alternative Hypothesis"),
    Conditional_Z_Score = pCond_Ridge,
    Conditional_aECS = pCond_aECS
  )
  return(res_raw)
}


experiment_configs <- list(
  list(scenario = "baseline", rho_seq = c(0.4, 0.9),  mismatch_prop = 0.00, noise_sd = 0.00,
       title_suffix = "",                     file_suffix = ""),
  list(scenario = "baseline", rho_seq = c(0.8, 0.95), mismatch_prop = 0.00, noise_sd = 0.00,
       title_suffix = "",                     file_suffix = ""),
  list(scenario = "mismatch", rho_seq = c(0.4, 0.9),  mismatch_prop = 0.05, noise_sd = 0.00,
       title_suffix = ", mismatch: 5%",        file_suffix = "_mismatch"),
  list(scenario = "mismatch", rho_seq = c(0.8, 0.95), mismatch_prop = 0.05, noise_sd = 0.00,
       title_suffix = ", mismatch: 5%",        file_suffix = "_mismatch"),
  list(scenario = "ld_noise", rho_seq = c(0.4, 0.9),  mismatch_prop = 0.00, noise_sd = 0.05,
       title_suffix = ", LD noise: SD=0.05",   file_suffix = "_ld_noise"),
  list(scenario = "ld_noise", rho_seq = c(0.8, 0.95), mismatch_prop = 0.00, noise_sd = 0.05,
       title_suffix = ", LD noise: SD=0.05",   file_suffix = "_ld_noise")
)

run_experiment <- function(config) {
  set.seed(2026)
  N_sim <- 100000
  nSnpA_seq <- c(50, 200)
  rho_seq <- config$rho_seq

  Sample_N <- 10000
  MISMATCH_PROP <- config$mismatch_prop
  NOISE_SD <- config$noise_sd

  pCut <- c(0.01, 5E-5, 5E-8)
  results_list <- list()

  for(r in rho_seq) {
    for(nA in nSnpA_seq) {
      cat(sprintf("Running: Rho=%.2f, nSnpA=%d (Sign Flipped=%.1f%%, Noise SD=%.2f)...\n",
                  r, nA, MISMATCH_PROP*100, NOISE_SD))
      df_null <- run_sim_raw_p(N_sim, nA, rho_true=r, h2_B=0.0, N=Sample_N,
                               sign_mismatch_prop=MISMATCH_PROP, ld_noise_sd=NOISE_SD)
      df_power <- run_sim_raw_p(N_sim, nA, rho_true=r, h2_B=0.005, N=Sample_N,
                                sign_mismatch_prop=MISMATCH_PROP, ld_noise_sd=NOISE_SD)
      results_list[[length(results_list) + 1]] <- df_null
      results_list[[length(results_list) + 1]] <- df_power
    }
  }

  df_all_raw <- bind_rows(results_list)

  summary_list <- lapply(pCut, function(th) {
    df_all_raw %>%
      group_by(nSnpA, Rho, Type) %>%
      summarise(Threshold = th, Z_Score_Rate = mean(Conditional_Z_Score < th, na.rm = TRUE),
                aECS_Rate = mean(Conditional_aECS < th, na.rm = TRUE), .groups = 'drop')
  })
  df_summary <- bind_rows(summary_list) %>% arrange(Threshold, Rho, Type, nSnpA)
  cat("\n=== Table: Type I Error and Power at Varying Thresholds ===\n")
  df_summary_print <- df_summary %>% mutate(Threshold = sprintf("%g", Threshold))
  print(as.data.frame(df_summary_print), row.names = FALSE)

  df_qq <- df_all_raw %>%
    pivot_longer(cols = c("Conditional_Z_Score", "Conditional_aECS"), names_to = "Method", values_to = "P_value") %>%
    mutate(P_value = ifelse(P_value < 1e-30, 1e-30, P_value)) %>%
    group_by(Rho, Type, Method, nSnpA) %>% arrange(P_value) %>%
    mutate(Expected = -log10(ppoints(n())), Observed = -log10(P_value), Rho_Label = sprintf("LD (ρ = %g)", Rho)) %>%
    ungroup() %>% mutate(nSnpA_factor = factor(nSnpA, levels = nSnpA_seq))

  plot_specs <- expand.grid(
    Rho = rho_seq,
    Type = c("Null Hypothesis", "Alternative Hypothesis"),
    stringsAsFactors = FALSE
  )

  dir.create("figures", showWarnings = FALSE, recursive = TRUE)

  for (i in seq_len(nrow(plot_specs))) {
    rho_i <- plot_specs$Rho[i]
    type_i <- plot_specs$Type[i]

    df_plot <- df_qq %>%
      filter(Rho == rho_i, Type == type_i)

    short_type <- ifelse(type_i == "Null Hypothesis", "Null", "Alternative")
    short_title <- sprintf("ρ = %g, %s%s", rho_i, short_type, config$title_suffix)
    rho_tag <- gsub("\\.", "", sprintf("%.2f", rho_i))
    type_tag <- ifelse(type_i == "Null Hypothesis", "null", "alternative")

    p_qq <- ggplot(
      df_plot,
      aes(x = Expected, y = Observed, color = nSnpA_factor, shape = Method)
    ) +
      geom_point(alpha = 0.6, size = 1.8) +
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "black",
        linewidth = 0.8
      ) +
      scale_color_brewer(palette = "Set1") +
      scale_shape_manual(
        values = c(
          "Conditional_aECS" = 16,
          "Conditional_Z_Score" = 17
        )
      ) +
      labs(
        title = short_title,
        x = expression(Expected~~-log[10](italic(P))),
        y = expression(Observed~~-log[10](italic(P)))
      ) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      )

    ggsave(
      filename = sprintf(
        "figures/calibratedLDError_rho%s_%s%s.png",
        rho_tag,
        type_tag,
        config$file_suffix
      ),
      plot = p_qq,
      width = 7,
      height = 6,
      dpi = 400
    )
  }

  invisible(df_summary)
}

summary_results <- vector("list", length(experiment_configs))
for (i in seq_along(experiment_configs)) {
  config <- experiment_configs[[i]]
  summary_results[[i]] <- run_experiment(config) %>%
    mutate(Scenario = config$scenario)
}


make_export_table <- function(summary_data) {
  type1_data <- summary_data %>%
    filter(Type == "Null Hypothesis") %>%
    select(Threshold, Rho, nSnpA, Z_Score_Rate, aECS_Rate) %>%
    rename(
      `Type I Error (Z_Score)` = Z_Score_Rate,
      `Type I Error (aECS)` = aECS_Rate
    )

  power_data <- summary_data %>%
    filter(Type == "Alternative Hypothesis") %>%
    select(Threshold, Rho, nSnpA, Z_Score_Rate, aECS_Rate) %>%
    rename(
      `Power (Z_Score)` = Z_Score_Rate,
      `Power (aECS)` = aECS_Rate
    )

  inner_join(
    type1_data,
    power_data,
    by = c("Threshold", "Rho", "nSnpA")
  ) %>%
    arrange(Threshold, Rho, nSnpA) %>%
    mutate(
      Threshold = case_when(
        Threshold == 5E-8 ~ "5E-8",
        Threshold == 5E-5 ~ "5E-5",
        Threshold == 0.01 ~ "0.01",
        TRUE ~ format(Threshold, scientific = TRUE)
      ),
      across(
        c(
          `Type I Error (Z_Score)`,
          `Type I Error (aECS)`,
          `Power (Z_Score)`,
          `Power (aECS)`
        ),
        ~ sprintf("%.4f", .x)
      )
    )
}

all_summaries <- bind_rows(summary_results)

baseline_table <- all_summaries %>%
  filter(Scenario == "baseline") %>%
  make_export_table()

write.csv(
  baseline_table,
  file = "figures/calibratedLDError_summary.csv",
  row.names = FALSE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)