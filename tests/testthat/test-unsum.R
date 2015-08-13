context("unsum gives expected output")

test_that("neast hasn't changed", {
  
  expect_equal(neast,
     structure(list(nodeNames = structure(c(1L, 2L, 6L, 3L, 5L, 4L),
     .Label = c("Boro Bank", "C&R", "NYB", "Priory", "RG", "Whitby"),
     class = "factor"), outSum....rowSums.neastTrue. = c(17.2, 
     16.8, 6.3, 2, 1.6, 2), inSum....colSums.neastTrue. = c(13.6, 
     16.2, 8.5, 0, 4.8, 2.8)), .Names = c("nodeNames", "outSum", "inSum"),
     row.names = c(NA, -6L), class = "data.frame"))
  
}) 