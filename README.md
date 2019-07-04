 [![Travis build status](https://travis-ci.com/baruuum/relabelKL.svg?token=k7R3D8yhYkrGz6yc4eQf&branch=master)](https://travis-ci.com/baruuum/relabelKL)

## C++ implementation of the relabeling method proposed by Stephens(2000)

While the relabeling algorithm proposed in Stephens (2000) has been shown to perform well in dealing with the label switching phenomenon in Bayesian finite mixture models, it has been pointed out that it is computationally expensive. Currently, the `label.switching::stephens` function implements Stepen's method. Yet, as it is written in `R`, it is rather slow, which makes it impractical to use on MCMC output of moderate size.

The `relabelKL` package is simply a `C++` implementation of Stephen's relabeling algorithm. It is implemented using the [Armadillo library](http://arma.sourceforge.net/) and sourced via the `RcppArmadillo` package. If available, the relabeling algorithm will activate OpenMP to run
parts of the code in parallel.

## How to use the package

Using the `devtools` package, `relabelKL` can be directly installed from github:

``` r
devtools::install_github(baruuum/relabelKL)
library(relabelKL)
```

Given `S`\*`N`\*`K` array of MCMC samples, call it `x`, where `N` is the number of units/individuals who belong to `K` latent classes/categories/extreme types, and where `S` is the number of posterior draws, calling the `relabelMCMC` function will relabel the output by Stephen's KL-algrorithm. There are only two options that can be specified: the numebr of maximum iterations to try (`maxit`) and whether intermediate results should be printed (`verbose`).

Running
``` r
res = relabelMCMC(x, maxit = 100, verbose = T)
```
will return a list with four elements: 

1.`permuted`: the same object as `x` but with the columns of each posterior draw relabeled
2. `perms`: a `S \times K` matrix which contains, for each draw, the optimal permutation. That is, `s`th row of `perms` shows how the `s`th posterior draw was permuted, so that applying the permutation in `perms` to the initial array, `x`, will result in the relabeled array
3. `iterations`: the number of iterations for which the algorithm was run
4. `status`: which is `0` if the algorithm has successfully converged and `1` otherwise

There are often other parameters in the model that depend on the labeling of the latent classes / extreme types. The `permuteMCMC` can be used in this situation. Suppose that the object `y` is an three-dimensional array, where the last dimension correspond to `S` posterior draws. After running `relabelMCMC` and obtaining the optimal permutations (or for any other relabeling algorithm that returns the permutation mapping), calling
```r
y.relabeled = permuteMCMC(y, perms = res$perms, what = "cols")
```
will permute the either the rows or the columns of each of the `S` sub-arrays of `y` according to `res$perms`. It is important to notice that the `permuteMCMC` function will assume that the last index of a three dimensional array represents the posterior draws, so that entering an array of dimensions, say,  `N`\*`B`\*`S` will lead to undefined behavior (assuming `S` stands for the draws). When a matrix is passed to the `permuteMCMC` function, it is assumed that it has dimensions `S`\*`K` and, thus, the function will always permute the columns. Lastly, sometimes we want to permute not only the rows or the columns but both simultaneously (which happens when the parameter of interest is a square matrix). If so, the `what = "both"` option can be used.

## Functions to use with `rstan` objects

The package comes with one dataset called `mmsbm`, which is a `stanfit` object obtained from running a mixed membership stochastic blockmodel on simulated data. The model has two parameters: `pi`, a matrix the mixed membership vectors, where each row indicates the probability of individual `i = 1,2,..,N`, belonging to type `k = 1, 2, ..., K`, and `theta` the so-called image matrix, which is of dimensions `K`\*`K` and reflects the association tendencies between the `K` pure types. 

The data can be loaded as
```r
data(mmsbm)
```
The `relabelKL` package provides a wraper function to extract posterior samples from a `stanfit` object and rearrange the parameter draws so that it can be directly passed to the `relabelMCMC` function. Calling
```r
pi.arr = extract_n_combine(mmsbm, par = "pi")
```
will create an array with the first dimension equal to the number of post-warmup draws times the number of chains run, and the other dimensions are identical to the dimensions specified in the `Stan` progam. For example, if `gamma` was specified as a `matrix[L,M]` object, the two last dimensions of the extracted object will be `L` and `M`m This object, then, can be passed to `relabelMCMC` or `permuteMCMC`. For example,
```r
rel = relabelMCMC(pi.arr, maxit = 50L, verbose = T)
rel.pi.arr = rel$permuted
theta.arr = extract_n_combine(mmsbm, par = "theta")
rel.theta.arr = permuteMCMC(theta.arr, rel$perms, "both")
```
can be used to relabel the posterior draws of both `pi` and `theta`. To monitor the convergence of the relabeled posterior draws, we have to transform the relabled arrays into a form that can be passed to the `rstan::monitor` function. This can be done by using the `to_stan_array` function:
```r
rstan::monitor(to_stan_array(rel.pi.arr))
rstan::monitor(to_stan_array(rel.theta.arr))
```
Similarly, traceplot might be inspected by using
```r
array_traceplot(to_stan_array(rel.pi.arr), par = "pi[1,2]")
array_traceplot(to_stan_array(rel.theta.arr), par = "theta")
```
The call will produce a `ggplot` object containing the traceplot of the single element `pi[1,2]` of the `pi` matrix; the secon call will produce a traceplot of all elements in the `theta` parameter.

## References

Stephens, M. 2000. "Dealing with label Switching in mixture models," *Journal of the Royal Statistical Society Series B*, 62, 795-809.
