test_that("add_ints works", {
  expect_equal(add_ints(1, 2), 3)
  expect_equal(add_ints(-1, 1), 0)
})
