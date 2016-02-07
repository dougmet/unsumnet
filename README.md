# `unsumnet` Reconstructing networks from aggregated data.

<a href="https://travis-ci.org/dougmet/unsumnet">
<img title="Build Status Images" src="https://travis-ci.org/dougmet/unsumnet.svg">
</a>

Welcome to the `unsumnet` project.

This project is an R interface to a C++ project to solve aggregated network data using simulated annealing. It is in heavy development but feel free to try it out now. If you find reproducible problems that have not already been raised then please create a ticket on the issues page.

To install from R you will need devtools installed:
```{r}
install.packages("devtools")
```

and then run

```{r}
library(devtools)
install_github(repo="unsumnet", username="dougmet", build_vignettes=TRUE)
```

There's a vignette already to help you get started. If you used the `build_vignettes=TRUE` argument then you can do

```{r}
vignette("Introduction", package="unsumnet")
```

A copy of this is in the wiki page [Introduction](https://github.com/dougmet/unsumnet/wiki/Introduction)

