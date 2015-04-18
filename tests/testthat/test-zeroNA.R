context("Zeros and NA handling")

test_that("Try some inputs with zeros and NAs", {
  # Simulate some input
  nBanks <- 20
  set.seed(10)
  sizes <- ceiling(rexp(nBanks, 1/20))
  constraints <- jitter(t(abs(floor(sapply(sizes, rnorm, n=2, sd=3)))))
  
  # Set some NAs
  na1 <- c(1,4,8,10)
  na2 <- 10:15
  constraints[na1, 1] <- NA
  constraints[na2, 2] <- NA
  # Note 10 has both NA
  
  # Set some zeroes
  zero1 <- c(5,6,17,18)
  zero2 <- c(3:8, 17)
  constraints[zero1, 1] <- 0
  constraints[zero2, 2] <- 0
  # Note 17 has both zero
  
  # Now some potentially confusing numbers
  constraints[1, 1] <- 1.022E-12
  
  expect_equal(length(which(constraints[,1]==0)), length(zero1))
  expect_true(all(which(constraints[,1]==0) == zero1))
  
})