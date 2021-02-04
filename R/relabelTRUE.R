#' Relabel MCMC output by minimizing KL-divergence to true assignment 
#' probabilities
#' 
#' Relabels the membership vectors of a mixed membership model or 
#' mixed membership stochastic blockmodel, by minimizing the KL-divergence
#' to a priori known true labels.
#'  
#' @param x an \code{S}\*\code{N}\*\code{K} array of MCMC samples containing
#'   the the assignment probabilities to latent classes/extreme types, where
#'   \code{N} is the number of units/individuals, \code{K} the number of
#'   latent classes, and \code{S} the number of posterior samples
#' @param x_true matrix of dimension \code{N}\*\code{K} which contains the 
#'   true assignment/mixed-membership probabilities
#' @param log_p if TRUE, treats elements in x as log-probabilities
#' @param renormalize if TRUE, renormalizes the rows of \code{x} and \code{x_true}
#' @param nthreads number of threads to use in parallel calculations
#' @param verbose if true, prints KL-divergence to true probabilities before
#'   and after relabeling
#' @return Returns a list of two elements: \code{permuted}, the relabeled
#'   array and, \code{perms}, the permutation pattern used for each sample
#' @export
relabelTRUE = function(x, x_true, log_p = TRUE, renormalize = FALSE, nthreads = 0L, verbose = TRUE) {
    
    if (!is.array(x))
        stop("x has to be an array")
    
    if (length(dim(x)) != 3L)
        stop("x has to be a three-dimensional array")
    
    if (!all(dim(x)[2:3] == dim(x_true)))
        stop("dimension mismatch between x and x_true")
    
    if (!log_p) {
        
        x = log(x)
        x_true = log(x_true)
        
    }
    
    if (sum(is.infinite(x) | x > 0.0) > 0)
        stop("log-probabilities in x have to be finite and less than zero")
    if (sum(is.infinite(x_true) | x_true > 0.0) > 0)
        stop("log-probabilities in x_true have to be finite and less than zero")
    
    
    if (
        !renormalize && !isTRUE(
            all.equal(
                apply(x, 1L, function(w) apply(w, 1L, lse)), 
                matrix(0.0, nrow = dim(x)[1L], ncol = dim(x)[2L]),
                check.attributes = FALSE
            )
        )
    )
        stop("probabilities in the rows of x don't sum to one")
    
    if (
        !isTRUE(
            all.equal(
                apply(x_true, 1L, lse),
                rep(0.0, nrow(x_true)),
                check.attributes = F
            )
        )
    )
        stop("rows in x_true do not sum to one")
    
    # relabel
    res = relabel_true_log(
        lphi        = aperm(x, c(2, 3, 1)), 
        lphi_true   = x_true, 
        renormalize = renormalize, 
        nthreads    = nthreads,
        verbose     = verbose
    )
    
    # change to original scale of input
    if (!log_p)
        res$permuted = exp(res$permuted)
    
    # change indexing of permutations to start from one
    res$perms = res$perms + 1L 
    # change ordering of dimensions 
    res$permuted = aperm(res$permuted, c(3, 1, 2))
    
    # add attributes of original array
    attr(res$permuted, "par") = attr(x, "par")
    attr(res$permuted, "org.att") = attr(x, "org.att")
    
    # return object
    res
    
}