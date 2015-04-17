context("Stop Checking")

test_that("Check stops fire for bad inputs", {
  x <- c(v1=10, v2=20, v3=30)
  expect_error(processInput(x))
})