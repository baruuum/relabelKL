#' Extract parameters into a long array from stanfit object
#'
#' Extracts MCMC samples stored in a stanfit objects and creats an array
#' which combines posterior samples across the MCMC chains *without*
#' permuting the indices, *excluding* the warmup phase, and returning
#' an object with the first dimension equal to the number of saved
#' posterior samples and the other dimensions equal to those specified
#' in the original Stan program (e.g., if parameter \code{x} was 
#' specified as \code{matrix[K,K]}, the resulting array will have dimensions
#' \code{S}\*\code{K}\*\code{K}, where \code{S} is the number of 
#' post-warmup posterior samples.)
#'
#' @param fit a stanfit object
#' @param par parameter to extract
#' @return the function returns an array of the extracted samples
#' @export
extract_n_combine = function(fit, par) {

    if (!requireNamespace("rstan", quietly = T))
        stop("package 'rstan' has to be installed to work with stanfit objects")

    if (!requireNamespace("abind", quietly = T))
        stop("please install the 'abind' package for this function to work")

    # extract pivot as array
    p = rstan::extract(fit, pars = par, permute = FALSE, inc_warmup = FALSE)

    # reformat post-warmup samples as array
    p.list = lapply(
        1:dim(p)[2],
        function (n) {
            tmp.array = p[, n, ]
            dim(tmp.array) = c(nrow(tmp.array), fit@sim$dims_oi[[par]])
            tmp.array
        })

    # combine arrays across chains
    res = abind::abind(p.list, along = 1L)
    attr(res, "par") = par
    attr(res, "org.att") = attributes(p)
    
    res
    
}

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

#' Traceplot for Stan-form arrays
#'
#' Plots traceplot for Stan-form arrays
#'
#' @param x an array created by calling \code{extract_n_combine}
#' @param par name of the parameter to plot
#' @return the function returns a \code{ggplot} object of the traceplot
#' @export
array_traceplot = function(x, par) {
    
    if(!requireNamespace("ggplot2", quietly = T))
        stop("install the 'ggplot2' package to create traceplots")
    if(!requireNamespace("reshape2", quietly = T))
        stop("'reshape2' package need to be installed to create traceplots")
    
    # original dims
    org.dim = dim(x)
    
    # get pars to plot
    p.pars = grep(
        gsub("(\\W)", "\\\\\\1", par) , 
        dimnames(x)$parameters, 
        value = TRUE)
    
    if (length(p.pars) == 0L)
        stop("plotting parameter (par) not found")
    
    # extract parameters
    x = x[, , p.pars]
    
    if (length(dim(x)) == 2L) {
        
        df = data.frame(value = c(x), 
                        chain = as.character(
                            rep(1:ncol(x), each = org.dim[1])
                        ),
                        iter = rep(1:org.dim[1], org.dim[2])
        )
                            
        ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, col = chain)) + 
        ggplot2::geom_line() + 
        ggplot2::theme_bw() +
        ggplot2::scale_color_viridis_d(
            option = "D", 
            begin = .1, 
            end = .8, 
            alpha = .6) + 
        ggplot2::ggtitle(par)
        
    } else {
    
        # reshape & creat data frame
        dim(x) = c(dim(x)[1] * dim(x)[2], dim(x)[3])
        df = data.frame(x)
        
        # add names and chain & iter information
        names(df) = p.pars
        df$chain = as.character(rep(1:org.dim[2], each = org.dim[1]))
        df$iter  = rep(1:org.dim[1], org.dim[2])
        
        # reshape to long format
        df = reshape2::melt(df, id.vars = c('iter', 'chain'))
        
        ggplot2::ggplot(df, ggplot2::aes(x = iter, y = value, col = chain)) + 
            ggplot2::geom_line() + 
            ggplot2::facet_wrap(~variable) + 
            ggplot2::theme_bw() +
            ggplot2::scale_color_viridis_d(
                option = "D", 
                begin = .1, 
                end = .8, 
                alpha = .6)
        
    }

}    
    