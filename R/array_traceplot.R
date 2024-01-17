#' Traceplot for Stan-form arrays
#'
#' Plots traceplot for Stan-form arrays
#'
#' @param x an array created by calling \code{extract_n_combine}
#' @param par name of the parameter(s) to plot entered as a sting. 
#'   This has to be either 
#'   \itemize{
#'     \item the name of a single parameter, e.g., \code{"theta[1,2]"}
#'     \item the name of an matrix/array/vector, e.g., \code{"theta"}
#'     \item or a vector of single parameters, e.g., \code{c("theta[1,1]", "theta[1,2]")}
#'   }
#' @return the function returns a \code{ggplot} object of the traceplot
#' @export
array_traceplot = function(x, par) {
    
    if(!requireNamespace("ggplot2", quietly = T))
        stop("install the 'ggplot2' package to create traceplots")
    if(!requireNamespace("reshape2", quietly = T))
        stop("'reshape2' package need to be installed to create traceplots")
    
    if (!is.character(par))
        stop("par has to be a character (string) object")
    
    # original dims
    org.dim = dim(x)
    
    # get pars to plot
    if (length(par) == 1L) {
        
        p.pars = grep(
            gsub("(\\W)", "\\\\\\1", par) , 
            dimnames(x)$parameters, 
            value = TRUE)
        
    } else if (length(par) > 1L) {
        
        p.pars = dimnames(x)$parameters[match(par, dimnames(x)$parameters)]
        
    } else {
        
        stop("par has length zero")
        
    }
    
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
            # ggplot2::scale_color_viridis_d(
            #     option = "D", 
            #     begin = .1, 
            #     end = .8, 
            #     alpha = .6) + 
            ggplot2::labs(x = "Iteration", y = "") +
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
            ggplot2::facet_wrap(~ variable, scales = "free") + 
            ggplot2::theme_bw() +
            # ggplot2::scale_color_viridis_d(
            #     option = "D", 
            #     begin = .1, 
            #     end = .8, 
            #     alpha = .6) +
            ggplot2::labs(x = "Iteration", y = "") + 
            ggplot2::theme(
                strip.text = ggplot2::element_text(hjust = .1),
                strip.background = ggplot2::element_blank()
            )
        
    }
    
}    
