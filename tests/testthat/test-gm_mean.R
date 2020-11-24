set.seed(420)
vector <- stats::rnorm(1000, mean = 50)

test_that("gm_mean works without NA values", {
  expect_equal(gm_mean(vector), 49.97892)
})

test_that("gm_mean works ignores NA values only when asked", {
  expect_equal(gm_mean(c(NA, vector), na.rm = TRUE), 49.97892)
  expect_warning(gm_mean(c(NA, vector), na.rm = FALSE))
})
