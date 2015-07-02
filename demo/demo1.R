## TO BE CONVERTED TO VIGNETTE

nodeNames <- c("Boro Bank", "Cook and Rea", "Bank of Whitby",
               "Priory Capital", "North Yorks Bank", "Redcar Group")

neastTrue <- rbind(c(0, 12.0, 1.2, 0, 1, 0),
               c(8.0, 0, 5.0, 0, 2.3, 1.5),
               c(4.1, 2.2, 1.0, 0, 0, 0),
               c(0.5, 0, 0, 0, 0.5, 0),
               c(0, 0, 0.3, 0, 0, 1.3),
               c(0, 2.0, 0, 0, 0, 0))

dimnames(neastTrue) <- list(nodeNames, nodeNames)

neastTrue
rowSums(neastTrue)
colSums(neastTrue)
neast <- data.frame(nodeNames, outSums=rowSums(neastTrue), inSums=colSums(neastTrue))

fit <- unsum(neast, 12, verbose=TRUE)

library(igraph)

plotUnsum <- function(aw) {
  g <- graph.adjacency(aw, weighted = TRUE)
  V(g)$size <- 50
  V(g)$color <- 'gold'
  if(!is.null(dimnames(aw)))
    V(g)$label <- dimnames(aw)[[1]]
  lo <- layout.circle(g)
  sw <- sum(aw)/50
  plot(g, edge.width=E(g)$weight/sw, edge.curved=is.mutual(g), layout=lo,
       margin=c(0,0,0,0))
}

plotUnsum(fit$AW)
plotUnsum(ATrue)


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


