cat("==== R / Gurobi test ====\n")
cat("R version: ", R.version.string, "\n")
cat("R binary:  ", Sys.which("R"), "\n")
cat("GUROBI_HOME:", Sys.getenv("GUROBI_HOME"), "\n")
cat("GRB_LICENSE_FILE:", Sys.getenv("GRB_LICENSE_FILE"), "\n\n")

# Load your package (this forces the Gurobi shared library to be loaded)
suppressPackageStartupMessages({
  library(SparsePortfolioSelection)
})

cat("SparsePortfolioSelection loaded OK\n")

# Small MIQP test (2 assets, k=1)
mu <- c(0.05, 0.10)
Sigma <- matrix(c(0.04, 0.0,
                  0.0,  0.09), nrow = 2)

res <- mve_miqp_search(
  mu = mu,
  sigma = Sigma,
  k = 1,
  time_limit = 5,
  threads = 1,
  verbose = FALSE
)

cat("\nMIQP status:", res$status, "\n")
cat("Weights:\n")
print(res$weights)

cat("\n==== SUCCESS: Gurobi is linked and working ====\n")
