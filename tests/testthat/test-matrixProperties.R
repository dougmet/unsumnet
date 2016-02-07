context('Calculating matrix properties')

# Don't want to set the seed for this test, use the time in case anyone
# set the seed upstream

timeSeed <- as.integer(Sys.time())

set.seed(timeSeed)

test_that('calc_max_edges correctly computes maximum number of edges', {
  for (i in 1:1000) {
    
    # All edges on except diagonal
    a <- 1-diag(nrow = 8)
    
    # How many rows and columns to wipe out
    c0 <- sample(0:7, 1, prob=c(0.3,0.1,0.1,0.1,0.1,0.1,0.05,0.05))
    r0 <- sample(0:7, 1, prob=c(0.3,0.1,0.1,0.1,0.1,0.1,0.05,0.05))
    
    # Wipe out some rows and columns
    if(c0) a[ , sample(1:8,c0)] <- 0
    if(r0) a[sample(1:8,r0) , ] <- 0
    
    expect_equal(calc_max_edges(rowSums(a), colSums(a)), sum(a>0),
                 info=paste0("Random seed was ", timeSeed, ", i=",i))
  }
})