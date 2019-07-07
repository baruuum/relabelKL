 [![Travis build status](https://travis-ci.com/baruuum/relabelKL.svg?token=k7R3D8yhYkrGz6yc4eQf&branch=master)](https://Travis-ci.com/baruuum/relabelKL)

## C++ implementation of the relabeling method proposed by Stephens(2000)

While the relabeling algorithm proposed in Stephens (2000) has been shown to perform well in dealing with the label switching phenomenon in Bayesian finite mixture models, it has been pointed out that it is computationally expensive. Currently, the `label.switching::stephens` function implements Stephens' method. Yet, as it is written in `R`, it is rather slow, which makes it impractical to use on MCMC output of moderate size.

The `relabelKL` package is simply a `C++` implementation of Stephen's relabeling algorithm. It is implemented using the [Armadillo library](http://arma.sourceforge.net/) and sourced via the `RcppArmadillo` package. If available, the relabeling algorithm will activate OpenMP to run
parts of the code in parallel.

## How to use the package

Using the `devtools` package, `relabelKL` can be directly installed from github:

``` r
devtools::install_github(baruuum/relabelKL)
library(relabelKL)
```

Given `S`\*`N`\*`K` array of MCMC samples, call it `x`, where `N` is the number of units/individuals who belong to `K` latent classes/categories/extreme types, and where `S` is the number of posterior draws, calling the `relabelMCMC` function will relabel the output by Stephen's KL-algorithm. There are only two options that can be specified: the number of maximum iterations to try (`maxit`) and whether intermediate results should be printed (`verbose`).

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

## Relabeling when true assignment probabilities are known

If the "true" labels of a stochastic blockmodel or latent class model are known in advance, assigning each individual to their true class is straightforward. Yet, there are situations in which we want to make the assignment probabilities of each posterior draw as close as possible to a set of fixed/true probabilities. These situations arise, for example, when bootstrapping finite mixture models or when we want to compare your posterior samples with elsewhere published results. Still other situations are those where we have calculated MAP or MLE estimates of the membership vectors or class assignments and want to use them as a pivots to relable the MCMC samples. In these cases, the `relabelTRUE` function might be used as follows:

```r
rel.true = relabelTRUE(x = x, x.true = x.true, verbose = T)
```

where `x.true` is a `N`\*`K` matrix of "true" assignment probabilities / mixed membership vectors. This function will try to relabel each posterior draw in `x` as close as possible to `x.true` in terms of KL-distances. 

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
will create an array with the first dimension equal to the number of post-warmup draws times the number of chains run, and the other dimensions are identical to the dimensions specified in the `Stan` program. For example, if `gamma` was specified as a `matrix[L,M]` object, the two last dimensions of the extracted object will be `L` and `M`m This object, then, can be passed to `relabelMCMC` or `permuteMCMC`. For example,
```r
rel = relabelMCMC(pi.arr, maxit = 50L, verbose = T)
rel.pi.arr = rel$permuted
theta.arr = extract_n_combine(mmsbm, par = "theta")
rel.theta.arr = permuteMCMC(theta.arr, rel$perms, "both")
```
can be used to relabel the posterior draws of both `pi` and `theta`. To monitor the convergence of the relabeled posterior draws, we have to transform the relabeled arrays into a form that can be passed to the `rstan::monitor` function. This can be done by using the `to_stan_array` function:
```r
rstan::monitor(to_stan_array(rel.pi.arr))
rstan::monitor(to_stan_array(rel.theta.arr))
```
Similarly, traceplot might be inspected by using
```r
array_traceplot(to_stan_array(rel.pi.arr), par = "pi[1,2]")
array_traceplot(to_stan_array(rel.theta.arr), par = "theta")
```
The call will produce a `ggplot` object containing the traceplot of the single element `pi[1,2]` of the `pi` matrix; the second call will produce a traceplot of all elements in the `theta` parameter.

## Comparison to `label.switching::stephens`

A simple comparison with the `stephens` function of the `label.switching` package might offer some insights into the efficiency gains. We first draw some random vectors that sum to one using the Dirichlet distribution.

```r
# set seed and dimensions
set.seed(123)
N = 30
K = 4
S = 1000

# generate random draws from Dirichlet
alpha = runif(1, 1, 5) * runif(K)
x = replicate(S, gtools::rdirichlet(N, alpha))

# reshape array to S * N * K
x = aperm(x, c(3, 1, 2))
```

Now, rigorous benchmarking will be super time-consuming; but as the difference in speed are quite large, we might use the `Sys.time()` function to get a sense of the efficiency gain. First, we use the `stephens` function of the `label.switching` package:

```r
# fit label.switching::stephens
start_t = Sys.time()
res1 = label.switching::stephens(x)
end_t = Sys.time()
print(end_t - start_t)
```
    
    Time difference of 3.256903 mins
    
Next, we use the `relabelMCMC` function:

```r
# fit relabelKL::relabelMCMC
start_t = Sys.time()
res2 = relabelKL::relabelMCMC(x, maxit = 100L, verbose = FALSE)
end_t = Sys.time()
print(end_t - start_t)
```

    Time difference of 7.258507 secs
    
Notice that it took **minutes** to relabel the draws with the `stephens` function but only **seconds** with the `relabelMCMC` function. Lastly, we make sure that the results agree as well.

```r
sum(res1$permutations != res2$perms)
```

    [1] 0
    
Other tests that were run to ensure that the algorithm works correctly can be found in the `/tests` directory of this github repo.


## Possible future extensions

1. Right now, the algorithm uses a *brute force* approach by comparing all permutations between the labels to minimize the KL-divergence. If the number of classes/extreme types gets large, this approach might become inefficient since the number of comparisons grows at the order of `K!`. The problem of minimizing the KL-divergence can be, however, formulated as an integer assignment problem for which better algorithms exists. Incorporating these might produce gains in efficiency.
2. Other relabeling algorithms might be rewritten in `C++` for faster implementation. The translation of algorithms contained in the `label.switching` package, for example, would be quite straightforward.

## References

Stephens, M. 2000. "Dealing with label Switching in mixture models," *Journal of the Royal Statistical Society Series B*, 62, 795-809.
