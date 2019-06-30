 [![Travis build status](https://travis-ci.com/baruuum/relabelKL.svg?branch=master)](https://travis-ci.com/baruuum/relabelKL)

## C++ Implementation of the Relabeling Method Proposed by Stephens(2000)

While the relabeling algorithm proposed in Stephens (2000) has been shown to perform well in dealing with the label switching phenomenon in Bayesian finite mixture models, it has been pointed out that it is computationally expensive. Currently, the `label.switching` package is the only `R` library that implements Stepen's method. Yet, as it's written in the `R` language, it is very slow, which makes it unpractical to use on MCMC output of moderate size.

The `relabelKL` package is a `C++` implementation of  Stephen's relabeling algorithm. It is implemented using the [Armadillo library](http://arma.sourceforge.net/) and sourced via the `RcppArmadillo` package.  

## References

Stephens, M. 2000. "Dealing with label Switching in mixture models," *Journal of the Royal Statistical Society Series B*, 62, 795-809.