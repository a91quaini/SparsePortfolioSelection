# Test that invalid input produces an error
test_that("load_data rejects invalid type", {
  expect_error(load_data(type = 2), "Invalid 'type'")
})

# ----- Test for type = 0 (Portfolios) -----
# Call load_data for portfolios
returns <- load_data(type = 0)
test_that("load_data(type = 0) returns a matrix with correct dimensions", {
  expect_true(is.matrix(returns))
  expect_equal(nrow(returns), 176)
  # Since each dummy portfolio contributes one column (after dropping the Date column),
  # the total number of columns should equal the number of portfolio objects.
  expect_equal(ncol(returns), 317)
})

# ----- Test for type = 1 (Individual Assets / CRSP) -----
# Call load_data for CRSP returns
returns <- load_data(type = 1)
test_that("load_data(type = 1) returns a matrix with correct dimensions", {
  expect_true(is.matrix(returns))
  expect_equal(nrow(returns), 176)
  # Removing the Date column leaves 5 columns.
  expect_equal(ncol(returns), 200)
})
