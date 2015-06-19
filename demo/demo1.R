## TO BE CONVERTED TO VIGNETTE

constraints <- matrix(rep(c(10,6,6,4,2,1),2),ncol=2)
fit <- unsum(constraints, 12)

library(igraph)


g <- graph.adjacency(fit$AW, weighted = TRUE)

plot(g, edge.width=E(g)$weight, edge.curved=is.mutual(g))
