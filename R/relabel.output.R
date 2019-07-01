#' Relabel MCMC output
#' 
#' Relabel MCMC output from finite mixture models using Stephens KL algorithm
#' 
#' @param x a \eqn{N \times K \times S} array of MCMC samples containing the
#'   the assignment probabilities to latent classes/extreme types, where
#'    \eqn{N} is the number of units/individuals, \eqn{K} the number of
#'    latent classes, and \eqn{S} the number of samples
#' @param maxit number of maximum iterations to use in the algorithm
#' @param verbose if true, prints the intermediate results
#' @return returns a list of two elements: \code{relabeled}, the relabeled
#'   array, and \code{perms} the permutation pattern used for each sample
#' @export
relabel = function(x, maxit = 100, verbose = TRUE) {
    
    if (!is.array(x))
        stop("x has to be an array")
    
    if (length(dim(x)) != 3L)
        stop("x has to be a three-dimensional array")
    
    if (sum(x < 0 | x > 1) > 0)
        stop("all elements in x have to lie between zero and one")
    
    if (
        !isTRUE(
            all.equal(
                apply(x, 3, rowSums), 
                matrix(1.0, nrow = dim(x)[1], ncol = dim(x)[3])
            )
        )
    )
        stop("rows in x do not sum to one")
    
    # relabel
    res = relabel_kl(x, maxit, verbose)
    # change indexing to start from one
    res$perms = res$perms + 1L 
    # return object
    return(res)
    
 
}

        
        