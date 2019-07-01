#' Relabel MCMC output by KL-algorithm
#' 
#' Relabel MCMC output from finite mixture models using Stephens KL algorithm
#' 
#' @param x an \eqn{N \times K \times S} array of MCMC samples containing the
#'   the assignment probabilities to latent classes/extreme types, where
#'    \eqn{N} is the number of units/individuals, \eqn{K} the number of
#'    latent classes, and \eqn{S} the number of samples
#' @param maxit number of maximum iterations to use in the algorithm
#' @param verbose if true, prints the intermediate results
#' @return a list of four elements: \code{relabeled}, the relabeled
#'   array; \code{perms} the permutation pattern used for each sample;
#'   \code{iterations} the number of iterations that were needed to permute
#'   the array; and \code{status} with 0 = converged, 1 = not converged.
#' @export
relabelMCMC = function(x, maxit = 100, verbose = TRUE) {
    
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
    
    # if maxit was reached, throw warning
    if (res$iterations == maxit) {
        
        warning("maximum number of iteration were reached")
        res$status = 1
        
    } else {
        
        res$staus = 0
        
    }
    
    # return object
    return(res)
    
 
}

#' Relabel MCMC output based on Permutation Matrix
#' 
#' Relabel MCMC output from finite mixture models using the permutations 
#' found by some algorithm
#' 
#' @param x an array of MCMC samples, where the last index has to be equal
#'   to the number of posterior samples, \code{S}, and the second-to-last
#'   index equal to the number of latent classes/extreme types, \code{K}.
#' @param perms a \code{S} times \code{K} matrix of permutations needed to
#'   transformed the raw MCMC output to the permuted output, e.g., the 
#'   \code{perms} element from running \code{relabelMCMC}.
#' @return \code{x} permuted according to \code{perms}
#' @details Given a permutation that resolves the label-switching, this 
#'   function applies the mapping to another object. If the object to 
#'   permute is a three-dimensional array of dimension \code{N} times
#'   \code{K} times \code{S}, the function will treat \code{S} as the 
#'   number of posterior draws and will permute the columns of the sub-array
#'   which is of dimensions \code{N} times \code{K}. If \code{x} is a 
#'   matrix, the function assumes that \code{nrow(x) == S} and 
#'   \code{ncol(x) == K}.
#' @export
permuteMCMC = function(x, perms) {
    
    if (!is.array(x))
        stop("x has to be an array (permuteMCMC)")
    
    n.dim = length(dim(x))
    S = nrow(perms)
    K = ncol(perms)
    
    if ((n.dim == 3) & ((dim(x)[2] != K) | (dim(x)[3] != S)))
        stop("size mismatch (permuteMCMC)")
    
    if ((n.dim == 2) & ((dim(x)[2] != K) | (dim(x)[1] != S)))
        stop("size mismatch (permuteMCMC)")
        
    
    if (n.dim == 3) {
        
        res = lapply(1:S, function(s) {
            x[,perms[s,],s]
        })
        
        simplify2array(res)
        
    } else if (n.dim == 2) {
        
        t(
            sapply(1:S, function(s) {
                x[s, perms[s,]]
            })
        )
        
    } else {

        stop("n.dim not 2 or 3 (permuteMCMC)")
        
    }
    
}
