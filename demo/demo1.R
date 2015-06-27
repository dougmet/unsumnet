## TO BE CONVERTED TO VIGNETTE

set.seed(11)
constraints <- matrix(rep(c(10,6,6,4,2,1),2),ncol=2)
fit <- unsum(constraints, 12)

library(igraph)


g <- graph.adjacency(fit$AW, weighted = TRUE)

plot(g, edge.width=E(g)$weight, edge.curved=is.mutual(g))


# This is a neat example where maxEntropy is rubbish.
set.seed(33)
a <- matrix(sample(0:1, 100, prob = c(0.85,0.15), replace=TRUE), ncol=10)
diag(a) <- 0
awTrue <- a * abs(rnorm(100, 40, 15))

rs <- rowSums(awTrue)
cs <- colSums(awTrue)

awME <- maxEntropy(rs,cs)
View(awTrue)

hist(awME, freq=FALSE)
hist(awTrue, freq=FALSE)

gTrue <- graph.adjacency(awTrue, weighted = TRUE)
plot(gTrue)

gME <- graph.adjacency(awME, weighted = TRUE)
plot(gME)


