#' Posterior samples from mixed membership stochastic blockmodel
#' 
#' A \code{stanfit} object obtained from fitting a mixed memebership stochastic
#' blockmodel to simulated data. 
#' 
#' @format A \code{stanfit} object:
#' \describe{
#'   \item{}{Model has two parameters: the membership vectors \code{pi} and the image matrix \code{theta}}
#'   \item{}{Two MCMC chains were run on simulated data using NUTS, each sampling 500 samples after 500 iterations of warmup}
#'   \item{}{The Stan-code can be accessed by typing \code{mmsbm@stanmodel}}
#' }
"mmsbm"