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

    if (!inherits(fit, "stanfit"))
        stop("fit has to be a stanfit object")
    
    if (!is.character(par))
        stop("par has to be a character object")
    
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
