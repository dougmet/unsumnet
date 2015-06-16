context('Calculating matrix properties')

test_that('calcMaxEdges correctly computes maximum number of edges', {
  for (i in 1:10000) {
    
    # All edges on except diagonal
    a <- 1-diag(nrow = 8)
    
    # How many rows and columns to wipe out
    c0 <- sample(0:7, 1, prob=c(0.3,0.1,0.1,0.1,0.1,0.1,0.05,0.05))
    r0 <- sample(0:7, 1, prob=c(0.3,0.1,0.1,0.1,0.1,0.1,0.05,0.05))
    
    # Wipe out some rows and columns
    if(c0) a[ , sample(1:8,c0)] <- 0
    if(r0) a[sample(1:8,r0) , ] <- 0
    
    expect_equal(calcMaxEdges(rowSums(a), colSums(a)), sum(a>0))
  }
})