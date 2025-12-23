# emp_analysis_managed_portfolios_daily.R
#
# Out-of-sample analysis for US managed portfolios (daily).
# Theory-consistent design: h = 1, roll by 1, moments with denominator T_in.

## ---- thread control: must be at the very top ------------------------------
Nn <- 180L
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

library(SparsePortfolioSelection)

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
PANEL_TYPE <- "US"
MISSINGS   <- "median"
N_ASSETS   <- 290L # total = 277
RNG_SEED   <- 12345

T_IN_GRID  <- c(360L, 480L, 600L)  # in-sample window lengths (days)
ADD_MKT    <- TRUE
ADD_FACTORS<- TRUE

K_MIN <- 3L
K_STEP<- 3L
K_CAP <- N_ASSETS - 1L

METHOD <- "lars"
REFIT  <- FALSE
PARALLEL <- TRUE
NORMALIZE <- TRUE

OUT_DIR <- file.path("inst", "empirics", "results", "managed_portfolios_daily")
FIG_DIR <- file.path("inst", "empirics", "figures", "managed_portfolios_daily")
sps_make_output_dirs(OUT_DIR, FIG_DIR)

set.seed(RNG_SEED)

# ---------------------------------------------------------------------
# Load + subset panel
# ---------------------------------------------------------------------
ld <- load_data(
  type = PANEL_TYPE, missing = MISSINGS, path = "data",
  frequency = "daily", add_mkt = ADD_MKT, add_factors = ADD_FACTORS,
  shuffle = TRUE
)

panel <- sps_subset_panel(ld$returns, rf = ld$rf, n_assets = N_ASSETS)
R_all  <- panel$returns
rf_vec <- panel$rf
N      <- panel$N

# if rf is missing but you still want turnover computed, set rf=0
if (is.null(rf_vec)) rf_vec <- rep(0, nrow(R_all))

k_grid <- sps_make_k_grid(N, k_min = K_MIN, k_step = K_STEP, k_cap = K_CAP)

panel_tag   <- PANEL_TYPE
factors_tag <- if (ADD_FACTORS) "ff3" else "nofactors"
mkt_tag     <- if (ADD_MKT) "mkt" else "nomkt"
stem_base   <- sps_make_stem_base(
  METHOD,
  refit = REFIT,
  panel_tag,
  factors_tag,
  mkt_tag,
  N,
  h = 1L,
  normalize_weights = NORMALIZE
)

# ---------------------------------------------------------------------
# Solver factory (closure pattern)
# ---------------------------------------------------------------------
solver_factory <- function(T_in) {
  force(T_in)
  function(mu, sigma, k) {
    res <- mve_lars_search(
      mu = mu, sigma = sigma, n_obs = T_in, k = k,
      ridge_epsilon   = 0.0,
      tol_nnl         = 1e-10,
      normalize_weights   = NORMALIZE,
      normalization_type  = 1L,
      use_refit       = REFIT,
      do_checks       = FALSE
    )
    list(weights = res$weights, selection = res$selection, status = res$status)
  }
}

# ---------------------------------------------------------------------
# Run suite
# ---------------------------------------------------------------------
suite <- run_empirical_suite_h1(
  R_excess = R_all, rf = rf_vec,
  T_in_grid = T_IN_GRID,
  k_grid = k_grid,
  solver_factory = solver_factory,
  out_dir = OUT_DIR,
  stem_base = stem_base,
  annualize = TRUE,
  return_paths = FALSE,
  parallel = PARALLEL,
  n_cores = Nn,
  parallel_backend = "fork"
)

# ---------------------------------------------------------------------
# Quick console tables
# ---------------------------------------------------------------------

# OOS Sharpe (annualized)
print_results(
  suite$k_grid,
  suite$mats$sharpe_oos,
  method_labels = suite$labels,
  digits = 4
)

# Cardinality check: for each k, how many windows had nnz(weights) != k
# (one column per T_in)
cat("\n=== nnz mismatch counts (per k, per T_in) ===\n")
print_results(
  suite$k_grid,
  suite$checks$nnz_mismatch_count,
  method_labels = suite$labels,
  digits = 1
)

# ---------------------------------------------------------------------
# Plots (+ PNG saving)
# ---------------------------------------------------------------------
plot_empirical_suite_h1(suite, fig_dir = FIG_DIR, stem_base = stem_base)
