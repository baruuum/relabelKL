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
        
        if (
            (p_dim == "cols" && dim(x)[3] != K) |
            (p_dim == "rows" && dim(x)[2] != K) | 
            (dim(x)[1] != S)
        )
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
