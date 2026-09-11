# =====================================================================
# simulation.R
# One Monte Carlo replication for
#   Jiang & Zhu (2026), "Extrapolating Treatment Effects in Multiple-
#   Score Regression Discontinuity Designs", 2026 AERA.
# Implements the design of Section 5: two CEF specifications M1, M2,
# 6 sample sizes, t in {0,1} controlling whether the CNIA holds.
# Output: model, t, n, p_ast, p_anova, est, sig.
# =====================================================================

simulation_fun <- function(sample_iter, model, t, n) {
  set.seed(sample_iter)
  
  ## 1. Potential-outcome CEFs ----------------------------------------
  if (model == 1) {                          # quadratic
    m0 <- function(s1, s2) 3*s1^2 + 3*s2^2 + t*s1*s2
    m1 <- function(s1, s2) 4*s1^2 + 4*s2^2 + t*s1*s2
  } else {                                   # Lee (2008) 5th-order
    P0 <- function(s) 1.27*s + 7.18*s^2 + 20.21*s^3 + 21.54*s^4 + 7.33*s^5
    P1 <- function(s) 0.84*s - 3.00*s^2 + 7.99*s^3 - 9.01*s^4 + 3.56*s^5
    I  <- function(s1, s2) t*(s1*s2 + 0.45*s1^2*s2 - 0.3*s1*s2^2)
    m0 <- function(s1, s2) 0.96 + P0(s1) + P0(s2) + I(s1, s2)
    m1 <- function(s1, s2) 1.04 + P1(s1) + P1(s2) + I(s1, s2)
  }
  
  ## 2. Data -----------------------------------------------------------
  x  <- runif(n, -0.25, 0.25)
  s1 <- 1.92*x + 0.52*(2*rbeta(n, 0.85, 1.15) - 1)     # U1
  s2 <- 1.92*x + 0.52*(2*rbeta(n, 0.85, 1.15) - 1)     # U2
  d  <- as.numeric(s1 >= 0 & s2 >= 0)                  # sharp and-case
  y  <- m0(s1, s2) + (m1(s1, s2) - m0(s1, s2))*d + rnorm(n, 0, 0.1295)
  
  ## 3. CNIA tests on the treated (Section 3.4) -----------------------
  f1 <- mgcv::gam(y ~ s(s1) + s(s2) + ti(s1, s2) + x,
                  subset = (d == 1), method = "REML")
  f0 <- mgcv::gam(y ~ s(s1) + s(s2) + x,
                  subset = (d == 1), method = "REML")
  p_ast   <- summary(f1)[[8]][3]                    # AST on ti(s1,s2)
  p_anova <- anova(f0, f1, test = "F")$`Pr(>F)`[2]  # ANOVA F-test
  
  ## 4. Boundary ATEs -> interior ATE at (0.5, 0.5), Eq. (3.8) --------
  rd <- try(rdmulti::rdms(Y = y, X = s1, C = c(0, 0.5, 0),
                          X2 = s2, C2 = c(0, 0, 0.5),
                          zvar = d, covs_mat = matrix(x, ncol = 1)),
            silent = TRUE)
  if (inherits(rd, "try-error")) {
    est <- sig <- NA_real_
  } else {
    est <- rd$Coefs[1, 3] + rd$Coefs[1, 2] - rd$Coefs[1, 1]
    sig <- sqrt(sum(rd$V[1, 1:3]))
  }
  
  c(model, t, n, p_ast, p_anova, est, sig)
}