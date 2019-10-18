 [![Travis build status](https://travis-ci.com/baruuum/relabelKL.svg?token=k7R3D8yhYkrGz6yc4eQf&branch=master)](https://Travis-ci.com/baruuum/relabelKL)

## C++ implementation of the relabeling method proposed by Stephens(2000)

While the relabeling algorithm proposed in Stephens (2000) has been shown to perform well in dealing with the label switching phenomenon in Bayesian finite mixture models and has the attractive feature of requiring no additional information beyond the posterior samples, it has been pointed out that it is computationally expensive. Currently, the `label.switching::stephens` function provides the only `R` implementation of Stephens' method. Yet, the function is rather slow, which makes it impractical to use on MCMC output of moderate size.

The `relabelKL` package is a `C++` implementation of Stephen's relabeling algorithm. It is implemented using the [Armadillo library](http://arma.sourceforge.net/) and sourced via the `RcppArmadillo` package. If available, the functions provided will use OpenMP to run parts of the code in parallel. The package provides also a set of helper functions to facilitate the integration with `stanfit` objects produced via the `rstan` package.

A comparison between the performance of the `label.switching::stephens` function and the `relabelKL::relabelMCMC` function is provided at the end of this document.

_Note_: On Mac OS, it sometimes happens that the installation fails with the error message that "-lgfortran" or "-lquadmath" cannot be found. Typing the following in your terminal will install the necessary libraries:
```
curl -OL http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```


## How to use the package

Using the `devtools` package, `relabelKL` can be directly installed from github:

``` r
devtools::install_github("baruuum/relabelKL")
library(relabelKL)
```

Relabeling MCMC output is done by the `relabelMCMC` function. The function takes as input an array of dimensions `S`\*`N`\*`K`, where `S` is the number of posterior draws, `N` the number of units/individuals, `K` the latent classes/categories/types to which they belong, and the entries are either probabilities or log-probabilities of unit `i = 1,2,...,N` belonging to one of the `k=1,2,..,K` classes/types. 

Three options might be specified: 

1. `maxit`: the number of maximum iterations (defaults to `maxit = 100`);
2. `verbose`: whether intermediate results should be printed (defaults to `verbose = TRUE`); and
3. `log.p`: whether the entries are log-probabilities or probabilities (defaults to `log.p = TRUE`)

Calling `relabelMCMC` on the array `x`
``` r
res = relabelMCMC(x)
```
returns a `list` of four elements: 

1. `permuted`: the same object as `x` but with the columns of each posterior draw relabeled
2. `perms`: a `S \times K` matrix which contains, for each draw, the optimal permutation. That is, the `s`th row of `perms` shows how the `s`th posterior draw was permuted, so that applying the permutation in `perms` to the initial array, `x`, will result in the relabeled array
3. `iterations`: the number of iterations for which the algorithm was run
4. `status`: which is `0` if the algorithm has successfully converged and `1` otherwise


## Relabeling other parameters of the model

There are often other parameters in the model that depend on the labeling of the latent classes / extreme types. The `permuteMCMC` function can be used in this situation. For an three-dimensional array, `y`, where the last dimension correspond to `S` posterior draws, the `permuteMCMC` can be called based on the results from the `relabelMCMC` output (or for any other relabeling algorithm that returns the permutation mapping). Calling, for example,
```r
y.relabeled = permuteMCMC(y, perms = res$perms, what = "dimension to permute")
```
will permute the either the rows, columns, or both the rows and columns of each of the `S` sub-arrays of `y` based on `res$perms`. 

1. Three strings might be passed to the `what` option: specifying `"rows`" will permute the rows of each `s=1,2,...,S` draw, specifying `"cols"` will permute the columns, and `"both"` will permute both the rows *and* columns.
1. It is important to notice that the `permuteMCMC` function assumes that the first index of a three dimensional array represents the posterior draws, so that entering an array of dimensions, say,  `N`\*`B`\*`S` will lead to undefined behavior. 
2. When a two-dimensional matrix, instead of a three-dimensional array, is passed to `permuteMCMC`, it is assumed that it has dimensions `S`\*`K` and, the function will always permute the columns. 

## Relabeling when true assignment probabilities are known

If the "true" labels of a stochastic blockmodel or latent class model are known in advance, assigning each individual to their true class is straightforward. Yet, there are situations in which we want to make the assignment probabilities of each posterior draw as close as possible to a set of known/true probabilities. These situations arise, for example, when MAP or MLE estimates of the membership vectors of mixed membership models are calculated and when we want to use them as a pivots to relable the MCMC samples. In these cases, the `relabelTRUE` function might be used as follows:

```r
rel.true = relabelTRUE(x = x, x.true = x.true, verbose = F, log.p = F)
```

where `x.true` is a `N`\*`K` matrix of "true" assignment probabilities / mixed membership vectors. This function will try to relabel each posterior draw in `x` as close as possible to `x.true` in terms of KL-distances. When the `log.p` option is set to `TRUE`, *both* `x` and `x.true` should be entered on the log-scale. Otherwise, the function will throw an error.

## Functions to use with `rstan` objects

The package comes with one dataset called `mmsbm`, which is a `stanfit` object obtained from running a mixed membership stochastic blockmodel on simulated data. The model has two parameters: `pi`, a matrix the mixed membership vectors, where each row indicates the probability of individual `i = 1,2,..,N` enacting type `k = 1, 2, ..., K` in the interaction process, and `theta` the so-called image matrix, which is of dimensions `K`\*`K` and reflects the association tendencies between the `K` pure types. 

The data can be loaded as
```r
data(mmsbm)
```
The `relabelKL` package provides a wraper function to extract posterior samples from a `stanfit` object and rearrange the parameter draws so that it can be directly passed to the `relabelMCMC` function. 

Calling
```r
pi.arr = extract_n_combine(mmsbm, par = "pi")
```
will combine post-warmup draws across all MCMC chains into the first dimension of a three-dimensional array, where the other dimensions are identical to the dimensions specified in the `Stan` program. For example, if `pi` was specified as a `matrix[L,M]` object, the first dimension of `pi.arr` will be of order `S`\*`J`, where `S` is as before the number of post-warmup draws and `J` the number of MCMC chains, and the two last dimensions of the extracted object will be `L` and `M`. This object, then, can be passed to `relabelMCMC` or `permuteMCMC`. 

Thus, by running
```r
rel = relabelMCMC(pi.arr, maxit = 50L, verbose = T, log.p = F)
rel.pi.arr = rel$permuted
```
`rel.pi.arr` will contain the relabled `pi` array and by running
```
theta.arr = extract_n_combine(mmsbm, par = "theta")
rel.theta.arr = permuteMCMC(theta.arr, rel$perms, "both")
```
we obtain the permuted image matrix. 

Lastly, to monitor the convergence of the relabeled posterior draws, we have to transform the relabeled arrays into a form that can be passed to the `rstan::monitor` function. This can be done by using the `to_stan_array` function:
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
res2 = relabelKL::relabelMCMC(x, maxit = 100L, verbose = FALSE, log.p = FALSE)
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

## References

Stephens, M. 2000. "Dealing with label Switching in mixture models," *Journal of the Royal Statistical Society Series B*, 62, 795-809.
