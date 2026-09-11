README
================================================================

Replication code for
"Extrapolating Treatment Effects in Multiple-Score Regression
Discontinuity Designs". 

----------------------------------------------------------------
FILES
----------------------------------------------------------------

simulation.R
  Defines simulation_fun(), which runs ONE Monte Carlo replication
  of the Section 5 design. Given (sample_iter, model, t, n), it
  generates the data, fits the two GAMs for the CNIA test, calls
  rdmulti::rdms() for the interior ATE, and returns
  (model, t, n, p_ast, p_anova, est, sig). Not meant to be run
  directly; it is sourced by run.R.

run_simulation.R
  Runs the full simulation grid: 10,000 replications x 2 models
  x 2 values of t x 6 sample sizes, in parallel. Saves the raw
  output to result_df.rds, then prints:
    - Table 1: empirical rejection rates of the AST and ANOVA tests
    - Table 2: RMSE, bias, CI length and coverage for the ATE
  Requires simulation.R in the same directory.

empirical.R
  Reproduces the Section 6 application on Colombia's Ser Pilo Paga
  program. Reads data/data_RD.dta, runs the CNIA / GCNIA tests,
  estimates the sharp ATT and fuzzy LATT, and draws the two
  heatmaps. Prints:
    - Table 3: p-values of the AST and ANOVA tests (sharp, fuzzy)
    - Table 4: ATT / LATT estimates with se, t-value, p-value
  and writes the txt grids and Figure2.pdf listed below.

----------------------------------------------------------------
OUTPUT FILES
----------------------------------------------------------------

output/Sharp_est.txt   11x11 grid of sharp ATE estimates.
output/Sharp_tvalue.txt
                       11x11 grid of corresponding t-values.
output/Fuzzy_est.txt   11x11 grid of fuzzy LATE estimates.
output/Fuzzy_tvalue.txt
                       11x11 grid of corresponding t-values.
output/Figure2.pdf     Figure 2: sharp and fuzzy heatmaps with a
                       shared legend.

----------------------------------------------------------------
USAGE
----------------------------------------------------------------

  source("run_simulation.R")        # simulation: Tables 1-2
  source("empirical.R")  # application: Tables 3-4, Figure 2

----------------------------------------------------------------
REQUIREMENTS
----------------------------------------------------------------

R >= 4.2 and the packages: mgcv, foreach, doParallel, xtable,
haven, rdmulti, lpdensity, ggplot2, dplyr, tidyr, RColorBrewer,
cowplot. rdmulti is not on CRAN; install it from the source listed
in the paper.

----------------------------------------------------------------
DATA
----------------------------------------------------------------

data_RD.dta comes from Londoño-Vélez, Rodríguez & Sánchez
(2020), AEJ: Economic Policy. Available at
https://www.aeaweb.org/journals/dataset?id=10.1257/pol.20180131