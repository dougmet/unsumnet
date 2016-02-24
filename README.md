<!-- README.md is generated from README.Rmd. Please edit that file -->
unsumnet
========

> Reconstructing networks from aggregated data.

[![Build Status](https://travis-ci.org/dougmet/unsumnet.svg?branch=master)](https://travis-ci.org/dougmet/unsumnet) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/dougmet/unsumnet?branch=master&svg=true)](https://ci.appveyor.com/project/dougmet/unsumnet) [![CRAN Status](http://www.r-pkg.org/badges/version/unsumnet)](http://www.r-pkg.org/pkg/unsumnet) [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/unsumnet)](http://www.r-pkg.org/pkg/unsumnet)

Welcome to the `unsumnet` project.

This project is an R interface to a C++ project to solve aggregated network data using simulated annealing. It is distributed as open source software under the MIT license. If you find reproducible problems that have not already been raised then please create a ticket on the issues page.

Installation
============

To install from R you will need devtools installed:

``` r
install.packages("devtools")
```

and then run

``` r
library(devtools)
install_github("dougmet/unsumnet", build_vignettes=TRUE)
```

There's a vignette already to help you get started. If you used the `build_vignettes=TRUE` argument then you can do

``` r
vignette("Introduction", package="unsumnet")
```

Introduction
============

There are many networks that we are interested in for which we only have aggregated data instead of the full structure. For example, in financial networks we often know the balance sheets (total assets and liabilities) of the individual institutions without knowing detailed exposure.

In matrix notation, the information we have is the row and column sums of the weighted adjacency network, \(AW_{ij}\).

\[ R_i = \sum_j AW_{ij} \] \[ C_i = \sum_j AW_{ji} \]

In general there are many possible networks that could satisfy these constraints. The aim of the `unsumnet` package is to use methods from statistical mechanics to evenly sample from these networks while allowing the user to choose the sparsity (number of edges).

Aggregated Data
---------------

Consider the following fictional banking system. It consists of six banks and they form a closed network of loans, which is represented by the weighted adjacency matrix, `neastTrue`,

``` r
library(unsumnet)
neastTrue
#>           Boro Bank  C&R Whitby NYB  RG Priory
#> Boro Bank       0.0 12.0    3.2   0 2.0    0.0
#> C&R             8.0  0.0    5.0   0 2.3    1.5
#> Whitby          4.1  2.2    0.0   0 0.0    0.0
#> NYB             1.5  0.0    0.0   0 0.5    0.0
#> RG              0.0  0.0    0.3   0 0.0    1.3
#> Priory          0.0  2.0    0.0   0 0.0    0.0
```

In reality this network is rarely known. It is much more common to have the aggregated data set, where the aggregation is the rowSums, `outSum=rowSums(neastTrue)` for outgoing loans and colSums `inSum=colSums(neastTrue)` for incoming loans.

``` r
neast
#>   nodeNames outSum inSum
#> 1 Boro Bank   17.2  13.6
#> 2       C&R   16.8  16.2
#> 3    Whitby    6.3   8.5
#> 4       NYB    2.0   0.0
#> 5        RG    1.6   4.8
#> 6    Priory    2.0   2.8
```

The primary function of this package is to find networks that fit these summed data. We wish to "unsum" the network.

Creating a network
------------------

The `unsum` function allows us to explore networks with the number of edges, `nEdges` as a control parameter. To make one network with 12 edges we simply call

``` r
set.seed(11)
fit <- unsum(neast, 12)
```

This returns an `unsumnet` object. The standard adjacency matrix is held in `fit$A` and the weighted adjacency matrix is `fit$AW`. The returned matrix will not be the same as the true network, but given the information we have it is equally likely because it satifies the constraints.

``` r
fit$AW
#>           Boro Bank   C&R Whitby NYB  RG Priory
#> Boro Bank       0.0 15.99    1.2   0 0.0   0.00
#> C&R             9.5  0.00    7.3   0 0.0   0.00
#> Whitby          1.6  0.00    0.0   0 2.8   1.86
#> NYB             1.1  0.00    0.0   0 0.0   0.94
#> RG              1.4  0.21    0.0   0 0.0   0.00
#> Priory          0.0  0.00    0.0   0 2.0   0.00
```

The `plot_unsum` method will plot a matrix and the S3 plot method will do the same for an `unsumnet` object.

``` r
plot_unsum(neastTrue, main="True network")
plot(fit, main="Reconstructed network")
```

<img src="inst/figs/README-unnamed-chunk-10-1.png" alt="True network compared to `unsumnet` reconstructed network"  /><img src="inst/figs/README-unnamed-chunk-10-2.png" alt="True network compared to `unsumnet` reconstructed network"  />
<p class="caption">
True network compared to `unsumnet` reconstructed network
</p>

Maximum Entropy Solution
------------------------

A commonly used solution for such aggregated data is known as the "maximum entropy" solution. The term is a little confusing from a statistical mechanics perspective but in this case it uses every possible edge and spreads out the weights as evenly as possible. The maximum entropy solution is provided by the `max_entropy` function. It takes the same inputs as `unsum`. In our example this looks quite different to the true network.

``` r
# Get the maximum entropy network
neastME <- max_entropy(neast)
# Compare to the true network
plot_unsum(neastTrue, main="True network")
plot_unsum(neastME, main="Max Entropy")
```

<img src="inst/figs/README-unnamed-chunk-11-1.png" alt="True network compared to maximum entropy network"  /><img src="inst/figs/README-unnamed-chunk-11-2.png" alt="True network compared to maximum entropy network"  />
<p class="caption">
True network compared to maximum entropy network
</p>

If no solution is possible then this function will return `FALSE`. This can be useful for testing impossible constraints.

Netted matrix
-------------

In stability terms it is often the netted positions that are more relevant. This is provided by the `netted_matrix` function. It keeps the positive elements \((AW-AW')_+\).

``` r
plot_unsum(netted_matrix(neastTrue))
```

![](inst/figs/README-unnamed-chunk-12-1.png)<!-- -->

Extra notes on `unsum`
----------------------

For large networks you might find the algorithm can take a very long time. The runtime is proportional to the number of edges. For very large networks with many edges it may become unusuable depending on the complexity of the constraints.

Another tricky case is when you ask for very few edges. This makes it harder to find solutions and the algorithm can become trapped in local minima (see below for how to prevent this). For example calling

``` r
unsum(neast, 9, verbose=TRUE) # won't converge
```

and will run around again and again.

The `unsum` function has many parameters with default values. These are detailed in the help page `?unsum` but the most common ones you might need to vary are below. `unsum` uses simulated annealing and this can take delicate tweaking of parameters if solutions are difficult to get.

-   `verbose`: Set this to `TRUE` if you want the algorithm to print to screen how its doing. Useful for big networks where it can take a long time to run.

-   `coolingRate`: This is an important parameter, it sets the rate of cooling. Each `mctSchedule` time steps the temperature is reduced by this factor. If it is too big then the quench will become trapped in a local minimum (an energy plateau). Too close to 1 and the simulation will take a long time.

-   `minError`: The mean squared error that the algorithm must get to before returning. This should be left alone unless you are really struggling to get a solution. If it's too big then the distribution risks becoming skewed.

-   `maxEdges`: If you set this to true *all* edges possible edges are switched on.

-   `noReturn`: If this is `TRUE` then no return edges are allowed. Essentially it's a netted matrix by construction.

Final note
----------

This package is package is under development and any feedback is welcome. If you find a reproducible problem then leave an issue on the [GitHub page](https://www.github.com/dougmet/unsumnet/). Other control parameters, such as controlling for known netted positions will be coming soon.

License
-------

MIT © [Douglas Ashton](https://github.com/dougmet)
