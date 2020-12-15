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
        !isTRUE(
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
#' @param what one of "rows", "cols," "both", or "matrix." Determines what dimensions
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
#'   and will permute the columns. When the \code{what} argument is set to "matrix," 
#'   it is assumed that the entered array is of dimension \code{S}\*\code{K}\*\code{L}\*\code{L}, and permutes the whole \code{L}\*\code{L} matrix across the \code{K} classes for each posterior draw. 
#' @export
permuteMCMC = function(x, perms, what = "cols") {
    
    p_dim = match.arg(
        what, 
        c("rows", "cols", "both", "matrix"),
        several.ok = FALSE
    )
    
    if (!is.array(x))
        stop("x has to be an array (permuteMCMC)")
    
    n_dim = length(dim(x))
    S = nrow(perms)
    K = ncol(perms)
    
    if (n_dim == 2) {
        
        if ((dim(x)[2] != K) | (dim(x)[1] != S))
            stop("size mismatch (permuteMCMC)")
        
        if (p_dim != "cols")
            stop("for two dimensional arrays, 'what' has to be set to 'cols'")

    }
    
    if (n_dim == 3) {
        
        if ((dim(x)[3] != K) | (dim(x)[1] != S))
            stop("size mismatch (permuteMCMC)")
        
        if (p_dim == "both" && dim(x)[2] != dim(x)[3]) 
            stop("parameter is not a square matrix but 'both' was specified")    
    }
    
    if (n_dim == 4) {
      
        if (p_dim != "matrix")
            stop("Only available option for the 'what' argument is 'matrix' for 4 dim arrays")
        
        stopifnot(
            dim(x)[2] == K,
            dim(x)[3] == dim(x)[4]
        )
            
    } 
        
    
    if (n_dim == 2) {
        
        res = t(sapply(1:S, function(s) { x[s, perms[s, ]] }))
        
    } else if (n_dim == 3) {
        
        res = switch(
            p_dim,
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
        
    } else if (n_dim == 4) {
        
        res = lapply(1:S, function(s) x[s, perms[s, ], ,])
        res = aperm(simplify2array(res), c(4, 1:3))
        
    } else {

        stop("n_dim must be smaller than 5 (permuteMCMC)")
        
    }
    
    # carry over attributes
    attr(res, "par") = attr(x, "par")
    attr(res, "org.att") = attr(x, "org.att")
    
    res
    
}
