if(!require(MASS)) install.packages("MASS")
if(!require(nnls)) install.packages("nnls")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(patchwork)) install.packages("patchwork")

library(MASS)
library(nnls)
library(ggplot2)
library(patchwork)

setwd("/public3/ly/SDESE_summary/aECS")
set.seed(2026)

p <- 30
rho <- 0.8
n_sim <- 100000

# 1. Structural decoupling and weight extraction

R <- rho^abs(outer(1:p, 1:p, "-"))
W <- coef(nnls(abs(R), rep(1, p)))
d_eff <- sum(W)

# 2. Extract the weighted kernel eigenvalue spectrum mu

W_sqrt <- diag(sqrt(W))
M <- W_sqrt %*% R %*% W_sqrt
mu <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
mu <- mu[mu > 1e-10]

# 3. Empirical sampling of the S_eff statistic

Z_matrix <- mvrnorm(n_sim, mu = rep(0, p), Sigma = R)
S_eff <- as.vector((Z_matrix^2) %*% W)

# 4. Exact spectral mixture benchmark sampling (Exact Spectral Mixture: sum mu_k * Chisq_1)

exact_spec_samples <- as.vector(matrix(rchisq(n_sim * length(mu), df = 1), nrow = n_sim) %*% mu)

# Figure A: Agreement between S_eff and the exact spectral mixture distribution (Q-Q Plot)

prob_seq <- seq(0.001, 0.999, length.out = 1000)
qq_data <- data.frame(
Theoretical = quantile(exact_spec_samples, probs = prob_seq),
Empirical   = quantile(S_eff, probs = prob_seq)
)

plotA <- ggplot(qq_data, aes(x = Theoretical, y = Empirical)) +
geom_point(color = "#1b9e77", alpha = 0.4, size = 1.5) +
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
theme_bw(base_size = 11) +
labs(title = "A: Exact Spectral Mixture Alignment",
subtitle = "Empirical S_eff vs Exact Spectrum (sum mu_k * Chisq_1)",
x = "Exact Theoretical Quantiles",
y = "Empirical S_eff Quantiles")

# Figure B: Type I Error control at stringent genome-wide significance levels

alpha_levels <- c(0.05, 0.01, 0.005, 0.001, 0.0005)

# Extract the theoretical rejection threshold at each alpha level from the exact null distribution

crit_thresholds <- quantile(exact_spec_samples, probs = 1 - alpha_levels)

# Calculate the empirical rejection rate from 100,000 samples under the null hypothesis

emp_rejections <- sapply(crit_thresholds, function(th) mean(S_eff > th))

type1_df <- data.frame(
Nominal_Log10   = -log10(alpha_levels),
Empirical_Log10 = -log10(emp_rejections)
)

plotB <- ggplot(type1_df, aes(x = Nominal_Log10, y = Empirical_Log10)) +
geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
geom_line(color = "#7570b3", linewidth = 1.2) +
geom_point(color = "#7570b3", size = 3) +
theme_bw(base_size = 11) +
labs(title = "B: Type I Error Control at Stringent Alpha",
subtitle = "Observed vs Expected rates down to alpha = 5 x 10^-4",
x = "-log10(Nominal Alpha)",
y = "-log10(Empirical Type I Error)")

# Figure C: Changes in the effective dimension d_eff and higher-order skewness calibration factor

rho_grid <- seq(0.1, 0.95, by = 0.05)
res_dynamics <- data.frame()

for(r in rho_grid) {
R_g <- r^abs(outer(1:p, 1:p, "-"))
W_g <- coef(nnls(abs(R_g), rep(1, p)))

M_g <- diag(sqrt(W_g)) %*% R_g %*% diag(sqrt(W_g))
mu_g <- eigen(M_g, symmetric = TRUE, only.values = TRUE)$values
mu_g <- mu_g[mu_g > 1e-10]

# Spectral skewness coefficient: skewness = sqrt(8) * sum(mu^3) / (sum(mu^2))^1.5

skew <- sqrt(8) * sum(mu_g^3) / (sum(mu_g^2)^1.5)

res_dynamics <- rbind(res_dynamics, data.frame(
rho = r,
d_eff = sum(W_g),
Skewness = skew
))
}

# Set the scaling factor for the secondary axis

scale_factor <- 5

plotC <- ggplot(res_dynamics, aes(x = rho)) +
geom_line(aes(y = d_eff, color = "Effective Dimension (d_eff)"), linewidth = 1.2) +
geom_point(aes(y = d_eff, color = "Effective Dimension (d_eff)"), size = 2) +
geom_line(aes(y = Skewness * scale_factor, color = "Spectral Skewness"), linewidth = 1.2, linetype = "dotdash") +
geom_point(aes(y = Skewness * scale_factor, color = "Spectral Skewness"), size = 2) +
scale_y_continuous(
name = "Effective Dimension (d_eff)",
sec.axis = sec_axis(~ . / scale_factor, name = "Spectral Skewness (gamma_1)")
) +
scale_color_manual(values = c("Effective Dimension (d_eff)" = "#d95f02",
"Spectral Skewness" = "#377eb8")) +
theme_bw(base_size = 11) +
labs(title = "C: Dimension Compression vs. Skewness Dynamics",
subtitle = "LD redundancy compression triggers higher spectral skewness",
x = "Linkage Disequilibrium (rho)",
color = "") +
theme(legend.position = c(0.68, 0.8),
legend.background = element_rect(fill = alpha("white", 0.6)))

# Arrange and export the minimalist 1x3 panel horizontally

final_three_plot <- plotA | plotB | plotC
ggsave("figures/S_eff.png", final_three_plot,
width = 16, height = 5, dpi = 800)
ggsave("figures/S_eff.pdf", final_three_plot,
width = 16, height = 5, device = grDevices::cairo_pdf, family = "Arial")