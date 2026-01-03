# emp_analysis_managed_portfolios_daily.R
#
# Out-of-sample analysis for US managed portfolios (daily).
# Theory-consistent design: h = 1, roll by 1, moments with denominator T_in.
# Uses: load_data(), sps_subset_panel(), sps_make_k_grid(), sps_make_stem_base(),
#       run_empirical_suite_h1(), print_results(), plot_empirical_suite_h1().

## ---- thread control: must be at the very top ------------------------------
N_CORES <- 180L
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

library(SparsePortfolioSelection)

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
DATA_TYPE <- "US"   # load_data() key
PANEL_TAG <- "US"   # for file stems/tags

MISSINGS <- "median"
N_ASSETS_REQUEST <- 290L   # total currently ~277; will cap automatically
RNG_SEED <- 12345

T_IN_GRID   <- c(360L, 480L, 600L)  # in-sample window lengths (days)
ADD_MKT     <- FALSE
ADD_FACTORS <- FALSE

K_MIN  <- 3L
K_STEP <- 3L
K_CAP_REQUEST <- NA_integer_        # NA => cap at N-1

METHOD    <- "lars"
REFIT     <- FALSE
NORMALIZE <- TRUE

ANNUALIZE <- TRUE
FREQUENCY <- "daily"

OUT_DIR <- file.path("inst", "empirics", "results", "managed_portfolios_daily")
FIG_DIR <- file.path("inst", "empirics", "figures", "managed_portfolios_daily")

.ensure_dir <- function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)
.ensure_dir(OUT_DIR)
.ensure_dir(FIG_DIR)

# Reproducibility under parallel
set.seed(RNG_SEED)
if (requireNamespace("parallel", quietly = TRUE)) {
  RNGkind("L'Ecuyer-CMRG")
}

# ---------------------------------------------------------------------
# Load + subset panel
# ---------------------------------------------------------------------
ld <- load_data(
  type = DATA_TYPE, missing = MISSINGS, path = "data",
  frequency = FREQUENCY, add_mkt = ADD_MKT, add_factors = ADD_FACTORS,
  shuffle = FALSE
)

if (is.null(ld$returns) || !is.matrix(ld$returns)) stop("load_data() did not return a returns matrix.")

R_raw   <- ld$returns
N_avail <- ncol(R_raw)

# Choose how many assets to use
if (is.na(N_ASSETS_REQUEST)) {
  N_ASSETS <- N_avail
} else {
  N_ASSETS <- min(as.integer(N_ASSETS_REQUEST), N_avail)
}
if (N_ASSETS < 2L) stop("Need at least 2 assets after subsetting.")

panel <- sps_subset_panel(R_raw, rf = ld$rf, n_assets = N_ASSETS)
R_all  <- panel$returns
rf_vec <- panel$rf
N      <- panel$N

# rf handling: if missing, set to 0 (still enables turnover computations)
if (is.null(rf_vec)) rf_vec <- rep(0, nrow(R_all))
if (length(rf_vec) != nrow(R_all)) stop("rf length must match number of rows in returns.")

# Explicit excess returns used for optimization/Sharpe
R_excess <- sweep(R_all, 1, rf_vec, "-")

# ---------------------------------------------------------------------
# k-grid + stems
# ---------------------------------------------------------------------
K_CAP <- if (is.na(K_CAP_REQUEST)) (N - 1L) else min(as.integer(K_CAP_REQUEST), N - 1L)
if (K_MIN > K_CAP) stop("K_MIN must be <= K_CAP.")

k_grid <- sps_make_k_grid(N, k_min = K_MIN, k_step = K_STEP, k_cap = K_CAP)

factors_tag <- if (ADD_FACTORS) "ff3" else "nofactors"
mkt_tag     <- if (ADD_MKT) "mkt" else "nomkt"

stem_base <- sps_make_stem_base(
  METHOD,
  refit = REFIT,
  panel_tag = PANEL_TAG,
  factors_tag = factors_tag,
  mkt_tag = mkt_tag,
  N = N,
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
      ridge_epsilon = 0.0,
      tol_nnl = 1e-10,
      normalize_weights = NORMALIZE,
      normalization_type = 1L,
      use_refit = REFIT,
      do_checks = FALSE
    )
    list(weights = res$weights, selection = res$selection, status = res$status)
  }
}

# ---------------------------------------------------------------------
# Parallel backend choice (portable)
# ---------------------------------------------------------------------
PAR_BACKEND <- "auto"

n_cores_eff <- N_CORES
if (requireNamespace("parallel", quietly = TRUE)) {
  n_cores_eff <- min(N_CORES, parallel::detectCores(logical = FALSE))
  n_cores_eff <- max(1L, n_cores_eff)
}

# ---------------------------------------------------------------------
# Run suite
# ---------------------------------------------------------------------
suite <- run_empirical_suite_h1(
  R_excess = R_excess,
  rf = rf_vec,
  T_in_grid = T_IN_GRID,
  k_grid = k_grid,
  solver_factory = solver_factory,
  out_dir = OUT_DIR,
  stem_base = stem_base,
  annualize = ANNUALIZE,
  frequency = FREQUENCY,
  return_paths = FALSE,
  parallel = TRUE,
  n_cores = n_cores_eff,
  parallel_backend = PAR_BACKEND,
  align_oos = "end"
)

# ---------------------------------------------------------------------
# Quick console tables
# ---------------------------------------------------------------------
cat("\n=== OOS Sharpe (", if (ANNUALIZE) "annualized" else "not annualized", ") ===\n", sep = "")
print_results(
  suite$k_grid,
  suite$mats$sharpe_oos,
  method_labels = suite$labels,
  digits = 4
)

cat("\n=== nnz mismatch counts (per k, per T_in) ===\n")
print_results(
  suite$k_grid,
  suite$checks$nnz_mismatch_count,
  method_labels = suite$labels,
  digits = 0
)

# ---------------------------------------------------------------------
# Save results (CSVs already saved by run_empirical_suite_h1)
# ---------------------------------------------------------------------
saveRDS(suite, file.path(OUT_DIR, paste0(stem_base, "_suite.rds")))

# ---------------------------------------------------------------------
# Plots (+ PNG saving)
# ---------------------------------------------------------------------
plot_empirical_suite_h1(suite, fig_dir = FIG_DIR, stem_base = stem_base)
