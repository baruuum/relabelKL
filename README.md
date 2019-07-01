 [![Travis build status](https://travis-ci.com/baruuum/relabelKL.svg?token=k7R3D8yhYkrGz6yc4eQf&branch=master)](https://travis-ci.com/baruuum/relabelKL)

## C++ Implementation of the Relabeling Method Proposed by Stephens(2000)

While the relabeling algorithm proposed in Stephens (2000) has been shown to perform well in dealing with the label switching phenomenon in Bayesian finite mixture models, it has been pointed out that it is computationally expensive. Currently, the `label.switching::stephens` function implements Stepen's method. Yet, as it is written in `R`, it is rather slow, which makes it impractical to use on MCMC output of moderate size.

The `relabelKL` package is simply a `C++` implementation of Stephen's relabeling algorithm. It is implemented using the [Armadillo library](http://arma.sourceforge.net/) and sourced via the `RcppArmadillo` package. So far, the provided function will **not** run in parallel, although the code could be parallelized in the future.

## How to use package

Using the `devtools` package, `relabelKL` can be directly installed from github:

``` r
devtools::install_github(baruuum/relabelKL)
library(relabelKL)
```

Given `N \times K\times S` array of MCMC samples, call it `x`, where `N` is the number of units/individuals who belong to `K` latent classes/categories/extreme types, and where `S` is the number of posterior draws, calling the `relabelMCMC` function will relabel the output by Stephen's KL-algrorithm. There are only two options that can be specified: the numebr of maximum iterations to try (`maxit`) and whether intermediate results should be printed (`verbose`).

Running
``` r
res = relabelMCMC(x, maxit = 100, verbose = T)
```
will return a list with four elements: 

1.`relabeled`: the same object as `x` but with the columns of each posterior draw relabeled;
2. `perms`: a `S\times K` matrix which contains, for each draw, the optimal permutation. That is, `s`th row of `perms` shows how the `s`th posterior draw was permuted, so that applying the permutation in `perms` to the initial array, `x`, will result in the relabeled array.
3. `iterations`: the number of iterations for which the algorithm was run
4. `status`: which is `0` if the algorithm has successfully converged and `1` otherwise.

There are often other parameters in the model that depend on the labeling of the latent classes / extreme types. The `permuteMCMC` can be used in this situation. After running `relabelMCMC` and obtaining the optimal permutations (or for any other relabeling algorithm that returns the permutation mapping), calling
```r
y.relabeled = permuteMCMC(y, res$perms)
```
will relabel `y` according to `res$perms`. It is important to notice that the `permuteMCMC` function will assume that a three dimensional array is of dimensions `N \times K \times S` and permutes the index corresponding to `K` (i.e, the columns of the "sub"-array of dimension `N \times K`), while it assumes that any matrix (or two-dimensional array) has dimensions `S \times K` and, thus, permutes the column  of that matrix.


## References

Stephens, M. 2000. "Dealing with label Switching in mixture models," *Journal of the Royal Statistical Society Series B*, 62, 795-809.