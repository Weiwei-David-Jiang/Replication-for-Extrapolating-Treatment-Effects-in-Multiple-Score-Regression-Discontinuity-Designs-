# =============================================================================
# empirical.R
# Empirical illustration for Jiang & Zhu (2026), Section 6.
# Colombia's Ser Pilo Paga (SPP) program.
# Produces Table 3 (CNIA / GCNIA tests), Table 4 (ATT / LATT) and Figure 2.
# =============================================================================

## ---- 0. Packages -----------------------------------------------------------
library(haven)
library(mgcv)
library(rdmulti)
library(lpdensity)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cowplot)

dir.create("output", showWarnings = FALSE)

## ---- 1. Data ---------------------------------------------------------------
df <- read_dta("data_RD.dta")

df <- df[!is.na(df$running_saber11), ]
df <- df[!is.na(df$running_sisben), ]
df <- df[, -which(names(df) %in% c(
  "running_saber11_placebo", "spadies_hq", "spadies_lq",
  "spadies_hq_pub", "spadies_lq_pub", "spadies_hq_pri", "spadies_lq_pri",
  "spadies_hq_id", "icfes_score_20132", "icfes_score_20142",
  "eligible_saber11", "sisben_score", "eligible_sisben", "icfes_per"))]

dat <- na.omit(df)
dat$r <- dat$running_saber11 / sd(dat$running_saber11)   # standardized SABER 11
dat$s <- dat$running_sisben  / sd(dat$running_sisben)    # standardized SISBEN
dat$d <- dat$beneficiary_spp
dat$y <- dat$spadies_any

# Restrict to [-3, 3] x [-3, 3] (paper: "within three standard deviations")
dat <- dat[dat$r <= 3 & dat$r >= -3 & dat$s <= 3 & dat$s >= -3, ]

# Column indices used throughout (see Section 5 of the paper for the design):
#   cols 4, 7:32       -> covariates passed to rdmulti::rdms()
#   cols 4, 7:34, 36   -> covariates + y, r, s used in the GAM tests
COVS_RDMS <- c(4, 7:32)
COVS_GAM  <- c(4, 7:34, 36)


# =============================================================================
# 2. Table 3 -- Testing CNIA (sharp) and GCNIA (fuzzy)
# =============================================================================
test_cnia <- function(dat, extra_filter = NULL) {
  idx <- dat$r <= 1 & dat$r >= -1 &
    dat$s <= 1 & dat$s >= -1 &
    dat$eligible_spp == 1
  if (!is.null(extra_filter)) idx <- idx & extra_filter
  
  dt <- dat[idx, COVS_GAM]
  preds <- setdiff(names(dt), c("y", "r", "s"))
  
  f_full <- gam(as.formula(paste("y ~ s(r) + s(s) + ti(r, s) +",
                                 paste(preds, collapse = " + "))),
                data = dt, method = "REML")
  f_red  <- gam(as.formula(paste("y ~ s(r) + s(s) +",
                                 paste(preds, collapse = " + "))),
                data = dt, method = "REML")
  anv <- anova(f_red, f_full, test = "F")
  
  c(AST   = round(summary(f_full)[[8]][3], 3),
    ANOVA = round(anv$`Pr(>F)`[2],       3))
}

tab3_sharp <- test_cnia(dat)
tab3_fuzzy <- test_cnia(dat, extra_filter = dat$beneficiary_spp == 1)

cat("Table 3  --  CNIA / GCNIA tests\n")
cat("  Sharp   AST =", tab3_sharp["AST"],
    ", ANOVA =", tab3_sharp["ANOVA"], "\n")
cat("  Fuzzy   AST =", tab3_fuzzy["AST"],
    ", ANOVA =", tab3_fuzzy["ANOVA"], "\n\n")


# =============================================================================
# 3. Figure 2 -- Heatmaps of interior ATE / LATE and t-values
# =============================================================================
knum <- 11
knot_point <- rbind(
  cbind(seq(0, 1, length.out = knum), 0),
  cbind(0, seq(0, 1, length.out = knum))[-1, ]
)

# ---- 3.1 Sharp scenario (ITT) ---------------------------------------------
fit_sharp <- rdms(Y = dat$y, X = dat$r,
                  C = knot_point[, 1], X2 = dat$s, C2 = knot_point[, 2],
                  zvar = dat$eligible_spp, covs_mat = dat[, COVS_RDMS])

# fit$Coefs[1:21], fit$V[1:21] are ordered as: boundary-1 (11 pts) then
# boundary-2 (10 pts, corner already counted). We assemble 11 x 11 grids.
est_mat <- se_mat <- matrix(NA_real_, knum, knum)
est_mat[1, ]      <- fit_sharp$Coefs[1:knum]              # tau(s1, 0)
est_mat[2:knum, 1] <- fit_sharp$Coefs[(knum + 1):(2 * knum - 1)]  # tau(0, s2)
se_mat[1, ]      <- sqrt(fit_sharp$V[1:knum])
se_mat[2:knum, 1] <- sqrt(fit_sharp$V[(knum + 1):(2 * knum - 1)])

for (i in 2:knum) for (j in 2:knum) {
  # Eq. (3.8): tau(s1,s2) = tau(s1,0) + tau(0,s2) - tau(0,0)
  est_mat[i, j] <- est_mat[i, 1] + est_mat[1, j] - est_mat[1, 1]
  # Eq. (3.10): Var = Var(s1,0) + Var(0,s2) + Var(0,0) for interior pts
  se_mat[i, j]  <- sqrt(se_mat[i, 1]^2 + se_mat[1, j]^2 + se_mat[1, 1]^2)
}
tv_mat_sharp <- est_mat / se_mat

write.table(est_mat, "output/Sharp_est.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(tv_mat_sharp, "output/Sharp_tvalue.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---- 3.2 Sharp ATT, Eq. (3.9) ---------------------------------------------
dat1 <- dat[dat$r >= 0 & dat$s >= 0 & dat$r <= 1 & dat$s <= 1, ]

fit_den_r <- lpdensity(data = dat1$r, grid = knot_point[, 1])
w_r <- fit_den_r$Estimate[2:knum, "f_p"]
t1  <- sum(w_r * fit_sharp$Coefs[2:knum]) / sum(w_r)
sig1 <- sum(w_r^2 * fit_sharp$V[2:knum]) / sum(w_r)^2

fit_den_s <- lpdensity(data = dat1$s, grid = knot_point[, 2])
w_s <- fit_den_s$Estimate[(knum + 1):(2 * knum - 1), "f_p"]
t2  <- sum(w_s * fit_sharp$Coefs[(knum + 1):(2 * knum - 1)]) / sum(w_s)
sig2 <- sum(w_s^2 * fit_sharp$V[(knum + 1):(2 * knum - 1)]) / sum(w_s)^2

att_sharp    <- t1 + t2 - fit_sharp$Coefs[1]
att_sharp_se <- sqrt(sig1 + sig2 + fit_sharp$V[1])
att_sharp_t  <- att_sharp / att_sharp_se
att_sharp_p  <- 2 * (1 - pnorm(abs(att_sharp_t)))

# ---- 3.3 Fuzzy scenario (per-protocol) -------------------------------------
fit_fuzzy <- rdms(Y = dat$y, X = dat$r,
                  C = knot_point[, 1], X2 = dat$s, C2 = knot_point[, 2],
                  zvar = dat$eligible_spp, covs_mat = dat[, COVS_RDMS],
                  fuzzy = dat$d)

est_mat <- se_mat <- matrix(NA_real_, knum, knum)
est_mat[1, ]      <- fit_fuzzy$Coefs[1:knum]
est_mat[2:knum, 1] <- fit_fuzzy$Coefs[(knum + 1):(2 * knum - 1)]
se_mat[1, ]      <- sqrt(fit_fuzzy$V[1:knum])
se_mat[2:knum, 1] <- sqrt(fit_fuzzy$V[(knum + 1):(2 * knum - 1)])

for (i in 2:knum) for (j in 2:knum) {
  est_mat[i, j] <- est_mat[i, 1] + est_mat[1, j] - est_mat[1, 1]
  se_mat[i, j]  <- sqrt(se_mat[i, 1]^2 + se_mat[1, j]^2 + se_mat[1, 1]^2)
}
tv_mat_fuzzy <- est_mat / se_mat

write.table(est_mat, "output/Fuzzy_est.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(tv_mat_fuzzy, "output/Fuzzy_tvalue.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---- 3.4 Fuzzy LATT, Eq. (4.4) ---------------------------------------------
# Weights omega(S1), omega(S2) estimated by linear regression of D on the
# running variable, over [0,1] x [0,1] where W = 1 and D = 1.
grid_r <- knot_point[2:knum, 1]   # 0.1, 0.2, ..., 1.0

lm_r <- lm(d ~ r, data = dat[dat$r >= 0 & dat$r <= 1 &
                               dat$s >= 0 & dat$s <= 1, ])
omega1 <- predict(lm_r, newdata = data.frame(r = grid_r)) /
  mean(dat$d[dat$r >= 0 & dat$r <= 1 & dat$s >= 0 & dat$s <= 1])

lm_s <- lm(d ~ s, data = dat[dat$r >= 0 & dat$r <= 1 &
                               dat$s >= 0 & dat$s <= 1, ])
omega2 <- predict(lm_s, newdata = data.frame(s = grid_r)) /
  mean(dat$d[dat$r >= 0 & dat$r <= 1 & dat$s >= 0 & dat$s <= 1])

fit_den_r <- lpdensity(data = dat1$r, grid = knot_point[, 1])
w_r <- fit_den_r$Estimate[2:knum, "f_p"]
t1  <- sum(w_r * omega1 * fit_fuzzy$Coefs[2:knum]) / sum(w_r)
sig1 <- sum((w_r * omega1)^2 * fit_fuzzy$V[2:knum]) / sum(w_r)^2

fit_den_s <- lpdensity(data = dat1$s, grid = knot_point[, 2])
w_s <- fit_den_s$Estimate[(knum + 1):(2 * knum - 1), "f_p"]
t2  <- sum(w_s * omega2 * fit_fuzzy$Coefs[(knum + 1):(2 * knum - 1)]) / sum(w_s)
sig2 <- sum((w_s * omega2)^2 * fit_fuzzy$V[(knum + 1):(2 * knum - 1)]) /
  sum(w_s)^2

latt_fuzzy    <- t1 + t2 - fit_fuzzy$Coefs[1]
latt_fuzzy_se <- sqrt(sig1 + sig2 + fit_fuzzy$V[1])
latt_fuzzy_t  <- latt_fuzzy / latt_fuzzy_se
latt_fuzzy_p  <- 2 * (1 - pnorm(abs(latt_fuzzy_t)))

cat("Table 4  --  ATT / LATT\n")
cat(sprintf("  Sharp   %.3f  (se %.3f,  t %.3f,  p %.3f)\n",
            att_sharp, att_sharp_se, att_sharp_t, att_sharp_p))
cat(sprintf("  Fuzzy   %.3f  (se %.3f,  t %.3f,  p %.3f)\n\n",
            latt_fuzzy, latt_fuzzy_se, latt_fuzzy_t, latt_fuzzy_p))


# =============================================================================
# 4. Figure 2 -- Two heatmaps + shared legend
# =============================================================================
make_heatmap_df <- function(est_file, tv_file) {
  est <- apply(read.table(est_file),   1:2, as.numeric)
  tv  <- apply(read.table(tv_file),    1:2, as.numeric)
  lab <- apply(est, 1:2, function(x) as.character(round(x, 3)))
  
  d_color <- as.data.frame(tv) |>
    tibble::rownames_to_column("y") |>
    pivot_longer(-y, names_to = "x", values_to = "color_value") |>
    mutate(x = as.numeric(sub("V", "", x)), y = as.numeric(y))
  d_label <- as.data.frame(lab) |>
    tibble::rownames_to_column("y") |>
    pivot_longer(-y, names_to = "x", values_to = "label") |>
    mutate(x = as.numeric(sub("V", "", x)), y = as.numeric(y))
  merge(d_color, d_label, by = c("x", "y"))
}

panel_theme <- theme_minimal() +
  theme(panel.grid       = element_blank(),
        axis.text.x      = element_text(size = 25, face = "italic",
                                        color = "black",
                                        margin = margin(t = 15), angle = 90),
        axis.text.y      = element_text(size = 25, face = "italic",
                                        color = "black",
                                        margin = margin(r = 20), vjust = 0.3),
        axis.title.y     = element_text(size = 40, angle = 90),
        axis.title.x     = element_text(size = 40, vjust = 2,
                                        margin = margin(t = 0.6, unit = "cm")),
        axis.ticks.length.x = unit(0.6, "cm"),
        axis.ticks.x     = element_line(color = "black", size = 3),
        axis.ticks.length.y = unit(0.6, "cm"),
        axis.ticks.y     = element_line(color = "black", size = 3))

make_panel <- function(df) {
  ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = color_value), color = "white", linewidth = 0.3) +
    geom_text(aes(label = label), color = "black", size = 5) +
    scale_fill_gradientn(colors = brewer.pal(5, "YlOrRd"),
                         limits = c(2, 28), breaks = c(2, 8, 15, 22, 28),
                         guide = "none") +
    scale_x_continuous(breaks = 1:11, labels = (0:10) * 0.1) +
    scale_y_continuous(breaks = 1:11, labels = (0:10) * 0.1) +
    coord_fixed() +
    labs(x = "Saber11", y = "SisBen") +
    panel_theme
}

df_sharp <- make_heatmap_df("output/Sharp_est.txt", "output/Sharp_tvalue.txt")
df_fuzzy <- make_heatmap_df("output/Fuzzy_est.txt", "output/Fuzzy_tvalue.txt")

p_sharp <- make_panel(df_sharp)
p_fuzzy <- make_panel(df_fuzzy)

legend_plot <- ggplot(df_sharp, aes(x, y, fill = color_value)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = brewer.pal(5, "YlOrRd"), name = "t-value",
    limits = c(2, 28), breaks = c(2, 8, 15, 22, 28),
    guide = guide_colorbar(
      barwidth  = unit(12, "cm"), barheight = unit(1, "cm"),
      title.position = "top",
      title.theme = element_text(size = 20),
      label.theme = element_text(size = 15, angle = 0, face = "italic"),
      ticks.linewidth = 1.5, ticks.colour = "grey20")) +
  theme(legend.position = "top", legend.justification = "center")
shared_legend <- get_legend(legend_plot)

pdf("output/Figure2.pdf", width = 16, height = 10)
plot_grid(shared_legend,
          plot_grid(p_sharp, p_fuzzy, nrow = 1, align = "h"),
          ncol = 1, rel_heights = c(0.15, 0.85))
dev.off()