
#' Collapse long-form array of posterior sample into short-form
#'
#' Collpases the long-form array created by \code{extract_n_combine} into
#' the array form that is created when calling \code{stan::extract} on a 
#' \code{stanfit} object with the options \code{permute = FALSE} and
#' \code{inc_warmup = FALSE}.
#'
#' @param x an array created by calling \code{extract_n_combine}
#' @return the function returns a reshaped array which is callable by
#'     \code{rstan::monitor} or \code{\link{array_traceplot}}
#' @export
to_stan_array = function(x) {
    
    if (!is.array(x))
        stop("x has to be an array")
    
    att.names = names(attributes(x))
    
    if (!("par" %in% att.names))
        stop("x needs parameter name ('par') as an attribute")
    
    if (!("org.att" %in% att.names))
        stop("x needs original attributes ('org.att') as an attribute")
    
    par = attr(x, "par")
    org.att = attr(x, "org.att")
    org.dims = org.att$dim
    
    res = array(NA, dim = org.dims)
    
    ix.s = 0L
    for (n in 1:org.dims[2]) {
        
        ix.s = ix.s + 1L
        ix.e = ix.s + org.dims[1] - 1L
        
        # reassign relabled pars
        if (length(dim(x)) == 2L) {
            
            tmp.arr = apply(x[ix.s:ix.e, ], 1L, c)
            
        } else if (length(dim(x)) == 3L) {
            
            tmp.arr = apply(x[ix.s:ix.e, , ], 1L, c)
            
        } else if (length(dim(x)) == 4L) {
            
            tmp.arr = apply(x[ix.s:ix.e, , ,], 1L, c)
            
        } else {
            
            stop("to_stan_array is currently not able to deal with arrays with more than 4 dimensions")
            
        }
        
        if (nrow(tmp.arr) != org.dims[3])
            stop("dimension mismatch")
        
        res[, n, ] = t(tmp.arr)
        
        ix.s = ix.e
        
    }
    
    dimnames(res) = org.att$dimnames
    res
    
}
