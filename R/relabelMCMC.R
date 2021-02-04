#' Relabel MCMC output by KL-algorithm
#' 
#' Relabel MCMC output from finite mixture models using Stephens KL algorithm
#' 
#' @param x an \code{S}\*\code{N}\*\code{K} array of MCMC samples containing
#'   the the assignment probabilities to latent classes/extreme types, where
#'   \code{N} is the number of units/individuals, \code{K} the number of
#'   latent classes, and \code{S} the number of posterior samples
#' @param log_p if TRUE, treats elements in x as log-probabilities
#' @param renormalize if TRUE, renormalizes the row of each slice of \code{x}
#' @param maxit number of maximum iterations to use in the algorithm
#' @param nthreads number of threads to use when relabeling (if OPENMP is enabled)
#' @param verbose if true, prints the intermediate results
#' @return Returns a list of four elements: \code{permuted}, the relabeled
#'   array; \code{perms} the permutation pattern used for each sample;
#'   \code{iterations} the number of iterations that were needed to permute
#'   the array; and \code{status} with 0 = converged, 1 = not converged.
#' @export
relabelMCMC = function(
    x,
    log_p = TRUE,
    renormalize = FALSE,
    maxit = 100, 
    nthreads = 0L,
    verbose = TRUE
) {
    
    if (!is.array(x))
        stop("x has to be an array")
    
    if (length(dim(x)) != 3L)
        stop("x has to be a three-dimensional array")
    
    if (!log_p) 
        x = log(x)
        
    if (sum(is.infinite(x) | x > 0.0) > 0)
        stop("log-probabilities have to be finite and less than zero")
    
    if (!renormalize && 
        !isTRUE(
            all.equal(
                apply(x, 1, function(w) apply(w, 1, lse)), 
                matrix(0.0, nrow = dim(x)[1], ncol = dim(x)[2]),
                check.attributes = FALSE
            )
        )
    )
        stop("probabilities in the rows of x don't sum to one")
    
    # relabel
    res = relabel_kl_log(
        lphi        = aperm(x, c(2, 3, 1)), 
        renormalize = FALSE,
        maxit       = maxit, 
        nthreads    = 0L,
        verbose     = verbose
    )
    
    # change to original scale of input
    if (!log_p)
        res$permuted = exp(res$permuted)
        
    # change indexing of permutations to start from one
    res$perms = res$perms + 1L 
    # change ordering of dimensions 
    res$permuted = aperm(res$permuted, c(3, 1, 2))
    
    # if maxit was reached, throw warning
    if (res$iterations == maxit) {
        
        warning("maximum number of iteration were reached")
        res$status = 1
        
    } else {
        
        res$status = 0
        
    }
    
    # add attributes of original array
    attr(res$permuted, "par") = attr(x, "par")
    attr(res$permuted, "org.att") = attr(x, "org.att")
    
    # return object
    res

}