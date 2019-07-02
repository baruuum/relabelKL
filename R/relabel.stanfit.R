#' Extract parameters into an array from stanfit object
#'
#' Extracts MCMC samples stored in a stanfit objects and creats an array
#' which combines posterior samples across the MCMC chains *without*
#' permuting the indices, *excluding* the warmup phase, and returning
#' an object with the first dimension equal to the number of saved
#' posterior samples and the other dimensions equal to those specified
#' in the Stan program
#'
#' @param fit a stanfit object
#' @param par parameter to extract
#' @param p.dims a named list which includes an element with the same name
#'    as \code{par} and which includes the dimensions of the parameter
#' @return an array of the extracted samples
extract_n_combine = function(fit, par, p.dims) {

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
            dim(tmp.array) = c(nrow(tmp.array), p.dims[[par]])
            tmp.array
        })

    # combine arrays across chains
    abind::abind(p.list, along = 1L)

}



#' Permute MCMC samples contained in a stanfit object using the KL-algorithm
#'
#' Permutes the samples stored within a stanfit object by relabeling them
#' using the KL-algorithm suggested by Stephens (2000)
#'
#' @param fit a stanfit object
#' @param par.list a list, each of which elements is a character vector of
#'     length 2. The first element of each entry of the list should be the
#'     name of the parameter and the second the dimensions that should be
#'     permuted (either "cols", "rows", or "both"). It is important that the
#'     first entry should specify the memebership/classification matrix
#'     based on which the KL-algorithm is run.
#' @param log.form a boolean that indicates whether the classification
#'     probabilities are estimated on the log-scale
#' @param maxit maximum number of iterations to be used in the KL-algorithm
#' @param verbose boolean that indicates whether intermediate messages
#'     should be printed out
#' @return a stanfit object with posterior samples permuted
#' @export
permute.stanfit = function(fit,
                       par.list,
                       log.form = FALSE,
                       maxit    = 100,
                       verbose  = FALSE) {

    if (!requireNamespace("rstan", quietly = T))
        stop("package rstan has to be loaded to work with stanfit objects")

    if (class(fit) != "stanfit")
        stop("fit has to be a stanfit object")

    if (fit@mode != 0)
        stop("fit@mode not equal zero")

    if (!is.list(par.list))
        stop("par.list has to be a list")

    if (sum(sapply(par.list, length) != 2L) > 0)
        stop("all elements of par.list have to be of length 2")

    if (sum(sapply(par.list, class) != "character") > 0)
        stop("all elements of par.list have to be characters")

    pivot = par.list[[1L]][1L]
    if (verbose)
        message(
            paste0("first element of par.list (",
                   pivot,
                   ") used as pivot to relabel MCMC ..."
            )
        )

    n.pars  = length(par.list)
    p.names = sapply(par.list, `[`, 1L)
    p.indx  = p.dims = vector("list", n.pars)
    names(p.indx) = names(p.dims) = p.names

    cc = 0
    for (pn in p.names) {

        cc = cc + 1

        if (!(pn %in% fit@sim$pars_oi))
            stop(paste0(pn, " was not sampled"))

        p.dims[[cc]] = fit@par_dims[[pn]]
        p.indx[[cc]] = which(
                grepl(
                    paste0('^', pn,'\\['),
                    fit@sim$fnames_oi
                    )
                )

        if (length(p.indx[[cc]]) == 0L)
            stop(paste0("couldn't find index for ", pn))

    }


    # get samp. params
    n.chains = length(fit@stan_args)
    n.warmup = fit@sim$warmup
    n.iter   = fit@sim$iter
    n.sample =  1 + (n.iter - n.warmup - 1) %/% fit@sim$thin

    if (verbose) message("extracting pivot ...")
    # extract pivot as array
    p.arr = extract_n_combine(fit, pivot, p.dims)

    if (verbose) message("relabeling pivot ...")
    # relabel combined array
    if (log.form) {

        res = relabelMCMC(apply(p.arr, 1:2, exp), maxit, verbose)

    } else {

        res = relabelMCMC(p.arr, maxit, verbose)

    }

    if (res$status != 0)
        stop(paste0("relabeling failed within ", maxit, " iterations"))

    # get permutations
    perms = res$perms

    # get permutated pivot
    permuted = res$permuted

    if (log.form) {

        permuted = apply(permuted, 1:2, log)

    }

    if (verbose) message("extracting and relabling rest of parameters ...")

    # extract other parameters and generate list
    perm.list = vector("list", n.pars)
    perm.list[[1L]] = permuted

    for (i in 2:n.pars) {

        perm.list[[i]] = permuteMCMC(
            extract_n_combine(fit, p.names[i], p.dims),
            perms,
            par.list[[i]][2L]
        )

    }

    if (verbose) message("reassigning MCMC to stanfit object ...")

    # reassign permuted pivot
    post.w = n.warmup + 1L
    ix.s = 0L
    for (n in 1:n.chains) {

        ix.s = ix.s + 1L
        ix.e = ix.s + n.sample - 1L

        # reassign relabled pars
        for (a in 1:n.pars) {

            tmp.arr = apply(perm.list[[a]][ix.s:ix.e, , ], 1L, c)

            if (nrow(tmp.arr) != length(p.indx[[a]]))
                stop(paste0("dimension mismatch, reassigning parameter (",
                            p.names[a], ")"))


            for (b in 1:length(p.indx[[a]])) {

                fit@sim$samples[[ n ]][[ p.indx[[a]][b] ]][post.w:n.iter] =
                    tmp.arr[b,]


            }

        }

        ix.s = ix.e

    }

    fit

}
