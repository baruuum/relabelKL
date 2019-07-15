#' Relabel MCMC output by KL-algorithm
#' 
#' Relabel MCMC output from finite mixture models using Stephens KL algorithm
#' 
#' @param x an \code{S}\*\code{N}\*\code{K} array of MCMC samples containing
#'   the the assignment probabilities to latent classes/extreme types, where
#'   \code{N} is the number of units/individuals, \code{K} the number of
#'   latent classes, and \code{S} the number of posterior samples
#' @param maxit number of maximum iterations to use in the algorithm
#' @param verbose if true, prints the intermediate results
#' @param log.p if TRUE, treats elements in x as log-probabilities.
#' @return Returns a list of four elements: \code{permuted}, the relabeled
#'   array; \code{perms} the permutation pattern used for each sample;
#'   \code{iterations} the number of iterations that were needed to permute
#'   the array; and \code{status} with 0 = converged, 1 = not converged.
#' @export
relabelMCMC = function(
    x, 
    maxit = 100, 
    verbose = TRUE,
    log.p = TRUE
) {
    
    if (!is.array(x))
        stop("x has to be an array")
    
    if (length(dim(x)) != 3L)
        stop("x has to be a three-dimensional array")
    
    if (log.p) {
        
        if (sum(is.infinite(x) | x > 0.0) > 0)
            stop("log-probabilities have to be finite and less than zero")
        
        if (
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
        res = relabel_kl_log(aperm(x, c(2, 3, 1)), maxit, verbose)
        
    } else {
        
        if (sum(x < 0 | x > 1) > 0)
            stop("all elements in x have to lie between zero and one")
        
        if (
            !isTRUE(
                all.equal(
                    apply(x, 1, rowSums), 
                    matrix(1.0, nrow = dim(x)[1], ncol = dim(x)[2]),
                    check.attributes = FALSE
                )
            )
        )
            stop("rows in x do not sum to one")
        
        # relabel
        res = relabel_kl(aperm(x, c(2, 3, 1)), maxit, verbose)
        
    }
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
#' @param x.true matrix of dimension \code{N}\*\code{K} which contains the 
#'   true assignment/mixed-membership probabilities
#' @param verbose if true, prints KL-divergence to true probabilities before
#'   and after relabeling
#' @return Returns a list of two elements: \code{permuted}, the relabeled
#'   array and, \code{perms}, the permutation pattern used for each sample
#' @export
relabelTRUE = function(x, x.true, verbose = FALSE) {
    
    if (!is.array(x))
        stop("x has to be an array")
    
    if (length(dim(x)) != 3L)
        stop("x has to be a three-dimensional array")
    
    if (sum(x < 0 | x > 1) > 0)
        stop("all elements in x have to lie between zero and one")
    
    if (sum(x.true < 0 | x.true > 1) > 0)
        stop("all elements in x.true have to lie between zero and one")
    
    if (
        !isTRUE(
            all.equal(
                apply(x, 1, rowSums), 
                matrix(1.0, nrow = dim(x)[1], ncol = dim(x)[2]),
                check.attributes = FALSE
            )
        )
    )
        stop("rows in x do not sum to one")
    
    if (
        !isTRUE(
            all.equal(
                rowSums(x.true),
                rep(1.0, nrow(x.true)),
                check.attributes = F
            )
        )
    )
        stop("rows in x.true do not sum to one")
    
    # relabel
    res = relabel_true(aperm(x, c(2, 3, 1)), x.true, verbose)
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



#' Relabel MCMC output based on Permutation Matrix
#' 
#' Relabel MCMC output from finite mixture models using the permutations 
#' found by some algorithm
#' 
#' @param x an array of MCMC samples, where the first index has to be equal
#'   to the number of posterior samples, \code{S}, and the last index equal 
#'   to the number of latent classes/extreme types, \code{K}.
#' @param perms a \code{S} times \code{K} matrix of permutations needed to
#'   transformed the raw MCMC output to the permuted output, e.g., the 
#'   \code{perms} element from running \code{relabelMCMC}.
#' @param what one of "rows", "cols," or "both". Determines what dimensions
#'   are permuted. For two-dimensional arrays with unequal dimensions
#'   (i.e., rectangular matrices), the columns are always permuted.
#' @return returns \code{x} permuted according to \code{perms}
#' @details Given a permutation that resolves the label-switching, this 
#'   function applies the mapping to another object. If the object to 
#'   permute is a three-dimensional array of dimension 
#'   \code{S}\*\code{N}\*\code{K}, the function will treat \code{S} as the 
#'   number of posterior draws and will permute either the columns or the
#'   rows of the sub-array which is of dimensions \code{N}*\code{K}. If
#'   \code{x} is a matrix, the function assumes that \code{nrow(x) == S} 
#'   and will permute the columns.
#' @export
permuteMCMC = function(x, perms, what) {
    
    p.dim = match.arg(what, 
                      c("rows", "cols", "both"),
                      several.ok = FALSE)
    
    if (!is.array(x))
        stop("x has to be an array (permuteMCMC)")
    
    n.dim = length(dim(x))
    S = nrow(perms)
    K = ncol(perms)
    
    if ((n.dim == 3) & ((dim(x)[3] != K) | (dim(x)[1] != S)))
        stop("size mismatch (permuteMCMC)")
    
    if ((n.dim == 2) & ((dim(x)[2] != K) | (dim(x)[1] != S)))
        stop("size mismatch (permuteMCMC)")
    
    if ((n.dim == 2) & p.dim %in% c("rows", "both")) 
        stop("if x is a two-dimensional array, make sure that the rows are the draws and the columns the labels")
    
    if ((n.dim == 3) & (p.dim == "both") & (dim(x)[2] != dim(x)[3])) 
        stop("parameter are not square matrices but 'both' was specified")
        
    
    if (n.dim == 3) {
        
        res = switch(
            p.dim,
            "rows" = lapply(1:S, function(s) {
                            x[s, perms[s, ], ]
                        }),
            "cols" = lapply(1:S, function(s) {
                            x[s, , perms[s, ]]
                        }),
            "both" = lapply(1:S, function(s) {
                            x[s, perms[s, ], perms[s, ]]
                        })
        )
        
        res = aperm(simplify2array(res), c(3, 1, 2))
        
    } else if (n.dim == 2) {
        
        res = t(sapply(1:S, function(s) { x[s, perms[s, ]] }))
        
    } else {

        stop("n.dim not 2 or 3 (permuteMCMC)")
        
    }
    
    # carry over attributes
    attr(res, "par") = attr(x, "par")
    attr(res, "org.att") = attr(x, "org.att")
    
    res
}
