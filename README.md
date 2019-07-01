 [![Travis build status](https://travis-ci.com/baruuum/relabelKL.svg?token=k7R3D8yhYkrGz6yc4eQf&branch=master)](https://travis-ci.com/baruuum/relabelKL)

## C++ Implementation of the Relabeling Method Proposed by Stephens(2000)

While the relabeling algorithm proposed in Stephens (2000) has been shown to perform well in dealing with the label switching phenomenon in Bayesian finite mixture models, it has been pointed out that it is computationally expensive. Currently, the `label.switching::stephens` function implements Stepen's method. Yet, as it is written in `R`, it is rather slow, which makes it impractical to use on MCMC output of moderate size.

The `relabelKL` package is simply a `C++` implementation of Stephen's relabeling algorithm. It is implemented using the [Armadillo library](http://arma.sourceforge.net/) and sourced via the `RcppArmadillo` package. So far, the provided function will **not** run in parallel, although the code could be parallelized in the future.

## References

Stephens, M. 2000. "Dealing with label Switching in mixture models," *Journal of the Royal Statistical Society Series B*, 62, 795-809.