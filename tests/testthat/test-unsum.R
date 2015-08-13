context("unsum gives expected output")

test_that("Expect that neast hasn't changed", {
  
  # If this is wrong then the system test will be wrong
  
  expect_equal(neast,
     structure(list(nodeNames = structure(c(1L, 2L, 6L, 3L, 5L, 4L),
     .Label = c("Boro Bank", "C&R", "NYB", "Priory", "RG", "Whitby"),
     class = "factor"), outSum....rowSums.neastTrue. = c(17.2, 
     16.8, 6.3, 2, 1.6, 2), inSum....colSums.neastTrue. = c(13.6, 
     16.2, 8.5, 0, 4.8, 2.8)), .Names = c("nodeNames", "outSum", "inSum"),
     row.names = c(NA, -6L), class = "data.frame"))
  
}) 

test_that("unsum gives expected output", {
  
  # Repeatable results
  set.seed(878)
  
  fit <- unsum(neast, nEdges=12, verbose=FALSE)
  
  # It's the right class
  expect_is(fit, "unsumnet")
  
  # Right number of elements
  expect_equal(length(fit),6)
  
  # Adjacency
  expect_equal(fit$A, structure(c(0L, 1L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L,
               1L, 1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L,
               0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L), .Dim = c(6L, 6L)))
  
  # Weighted adjacency
  expect_equal(fit$AW, structure(c(0, 6.26476847850537, 6.30000000380884, 0,
               1.03523151768579, 0, 14.2000000017013, 0, 0, 0, 0, 1.9999999982987,
               2.15135016872519, 5.78388134824892, 0, 0, 0.564768483025897, 0, 0,
               0, 0, 0, 0, 0, 0.848649818030312, 3.95135018196969, 0, 0, 0, 0,
               0, 0.799999997797397, 0, 2.0000000022026, 0, 0), .Dim = c(6L, 6L),
               .Dimnames = list(c("Boro Bank", "C&R", "Whitby", "NYB", "RG", "Priory"),
                                c("Boro Bank", "C&R", "Whitby", "NYB", "RG", "Priory"))))
  
  # Weights
  expect_equal(fit$W,
  structure(c(1e-05, 6.38149495443137, 6.47278201898688, 8.22236741140327, 
  1.09372318124607, 13.205127253305, 14.0561337193264, 1e-05, 0.27769050569087, 
  8.14702954153692, 16.1077949332101, 1.97084909469106, 2.10206435017032, 
  5.72852591616492, 1e-05, 19.0304050640966, 0.580158269290395, 
  8.69860716946537, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 0.84418385453581, 
  3.98421388320749, 12.8233974556312, 6.95533003265631, 1e-05, 
  6.8033043052448, 16.9131108877834, 0.678984371348396, 6.36157880490459, 
  2.00893103330332, 12.6803754397867, 1e-05), .Dim = c(6L, 6L)))
  
  # Reports
  expect_equal(fit$targetEdges, 12)
  
  expect_equal(fit$nEdges, 12)
  
  expect_equal(fit$Results, structure(c(1L, 0L, 0L, 0L),
                  .Names = c("Success", "TimeOut", "Plateau", "Invalid")))
  
  
})