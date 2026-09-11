# =====================================================================
# run_simulation.R
# Runs the full Monte Carlo grid in parallel and prints Table 1
# (rejection rates) and Table 2 (ATE at (0.5, 0.5)).
# =====================================================================
source("simulation.R")
library(foreach); library(doParallel); library(xtable)

## ---- 1. Parameter grid ---------------------------------------------
grid <- expand.grid(iter = 1:1e4, model = 1:2, t = 0:1,
                    n = c(300, 500, 1000, 1500, 2000, 3000))

## ---- 2. Parallel run ------------------------------------------------
cl <- makeCluster(max(parallel::detectCores() - 1L, 1L))
registerDoParallel(cl)
t0 <- Sys.time()

raw <- foreach(i = seq_len(nrow(grid)),
               .packages = c("mgcv", "rdmulti"),
               .export   = "simulation_fun") %dopar% {
                 tryCatch(with(grid[i, ], simulation_fun(iter, model, t, n)),
                          error = function(e) rep(NA_real_, 7))
               }
stopCluster(cl)
cat("Elapsed:", format(Sys.time() - t0), "\n")

result_df <- setNames(as.data.frame(do.call(rbind, raw)),
                      c("model", "t", "n", "test_ast", "test_anv",
                        "est", "sig"))
saveRDS(result_df, "result_df.rds")

## ---- 3. Table 1: rejection rates ------------------------------------
Ns <- c(300, 500, 1000, 1500, 2000, 3000)

rej <- function(nn, m, tt) {
  idx <- with(result_df, n == nn & model == m & t == tt)
  c(mean(result_df$test_ast[idx] < .05, na.rm = TRUE),
    mean(result_df$test_anv[idx] < .05, na.rm = TRUE))
}
tab1 <- t(sapply(Ns, function(nn)
  c(n = nn, rej(nn,1,0), rej(nn,1,1), rej(nn,2,0), rej(nn,2,1))))
colnames(tab1) <- c("n",
                    paste0("M1_t0_", c("AST","ANOVA")), paste0("M1_t1_", c("AST","ANOVA")),
                    paste0("M2_t0_", c("AST","ANOVA")), paste0("M2_t1_", c("AST","ANOVA")))
print(xtable(tab1, digits = 3))

## ---- 4. Table 2: ATE at (0.5, 0.5), t = 0 ---------------------------
stats <- function(idx, truth) {
  est <- result_df$est[idx]; sig <- result_df$sig[idx]
  c(RMSE = sqrt(mean((est - truth)^2, na.rm = TRUE)),
    Bias = mean(est - truth,           na.rm = TRUE),
    CL   = mean(2*1.96*sig,            na.rm = TRUE),
    CR   = mean(abs(est - truth) <= 1.96*sig, na.rm = TRUE))
}
tab2 <- t(sapply(Ns, function(nn)
  c(n = nn,
    stats(with(result_df, n == nn & model == 1 & t == 0), 0),
    stats(with(result_df, n == nn & model == 2 & t = 0), 0.08))))