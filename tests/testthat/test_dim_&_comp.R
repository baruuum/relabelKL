context("Check Dimensions and Compare with label.swithching package")

has_ls = requireNamespace("label.switching", quietly = T)

# function to draw from dirichlet distribution
rdirichlet = function(N, alpha, log.p = F) {
    
    l = length(alpha)
    x = matrix(rgamma(l * N, alpha), ncol = l, byrow = T)
    if (log.p) {
      log(x) - log(as.vector(x %*% rep(1, l)))
    } else {
      x/as.vector(x %*% rep(1, l))
    }
}

test_that("relabelMCMC throws appropriate errors with log/non-log input", {
  
    N = sample(20:100, 1)
    K = sample(2:4, 1)
    pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
    S = 100

    # generate array 
    test_ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    # reshuffle dimensions to match arg
    test_ar = aperm(test_ar, c(3, 1, 2))
    
    expect_error(relabelMCMC(test_ar, 100, F, T))
    expect_error(relabelMCMC(log(test_ar), 100, F, F))
  
})


test_that("relabelMCMC works in both verbose and non-verbose mode", {
  
    expect_output(relabelMCMC(test_ar, 100, T, F), regexp = ".")
    expect_silent(relabelMCMC(test_ar, 100, F, F))
    expect_output(relabelMCMC(log(test_ar), 100, T, T), regexp = ".")
    expect_silent(relabelMCMC(log(test_ar), 100, F, T))
    
})

test_that("results coincide with stephens function from the label.switching package", {
    
    # how many tests to run
    n_test = 5
    
    for (tt in 1:n_test) {
        
        # params
        N = sample(20:100, 1)
        K = sample(2:4, 1)
        pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
        S = 100

        # generate array 
        test_ar = simplify2array(
                    lapply(1:S, function(w) rdirichlet(N, pvec))
                  )
        # reshuffle dimensions to S * N * K
        test_ar = aperm(test_ar, c(3, 1, 2))
        # randomly relable some of the draws
        rel_draws = sample.int(S, floor(S/4), replace = F)
        for (s in rel_draws) 
            test_ar[s, , ] = test_ar[s, , sample.int(K, K, F)]

        # results from the label.switching package
        if (has_ls) 
            res1 = label.switching::stephens(test_ar)
        
        # results of package
        res2 = relabelMCMC(test_ar, maxit = 100, verbose = F, log.p = F)
        
        # check dimensions
        expect_true(identical(dim(res2$permuted), dim(test_ar)))
        # each label appears only once in each row
        expect_equal(
            sum(sapply(1:S, function(s) {
                length(res2$perms[s,]) != length(unique(res2$perms[s,]))
            })), 
            0
        )
        
        # compare with label.switching::stephens
        if (has_ls)
            expect_true(sum(res1$permutations!=res2$perms) == 0)
        
    }
    
})

test_that("relabeled array matches permutations", {
    
    # params
    N = sample(20:100, 1)
    K = sample(2:4, 1)
    pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
    S = 100

    # generate array 
    test_ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    # reshuffle dimensions to S * N * K
    test_ar = aperm(test_ar, c(3, 1, 2))
    # randomly relabel some of the draws
    rel_draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel_draws) 
        test_ar[s, , ] = test_ar[s, , sample.int(K, K, F)]
    
    expect_error(
        relabelMCMC(aperm(test_ar, c(2,3,1)), verbose = F, log.p = F)
    )

    res = relabelMCMC(test_ar, verbose = F, log.p = F)    
    
    expect_equal(permuteMCMC(test_ar, res$perms, "cols"), res$permuted)
    expect_error(permuteMCMC(test_ar, res$perms, "both"))
    
    samp_ind = sample.int(dim(test_ar)[2], 1L)
    expect_equal(
        permuteMCMC(test_ar[, samp_ind, ], res$perms, "cols"),
        res$permuted[, samp_ind,]
    )
    expect_error(permuteMCMC(test_ar[, samp_ind, ], res$perms, "rows"))
    expect_error(permuteMCMC(test_ar[, samp_ind, ], res$perms, "both"))
    
})

test_that("2nd relabeling results in identity mapping", {
    
    # params
    N = sample(20:100, 1)
    K = sample(2:4, 1)
    pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
    S = 100

    # generate array 
    test_ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    # reshuffle dimensions to S * N * K
    test_ar = aperm(test_ar, c(3, 1, 2))
    # randomly relabel some of the draws
    rel_draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel_draws) 
        test_ar[s, , ] = test_ar[s, , sample.int(K, K, F)]

    res = relabelMCMC(test_ar, maxit = 100, verbose = F, log.p = F)    
    res2 = relabelMCMC(res$permuted, maxit = 100, verbose = F, log.p = F)
    
    expect_equal(res2$iterations, 0L)
    expect_equal(res2$status, 0L)
    expect_equal(res2$permuted, res$permuted)
    expect_equal(res2$perms, matrix(rep(1:K, S), nr = S, byrow = TRUE))
    
    # same test for log-scale
    res_log = relabelMCMC(log(test_ar), maxit = 100, verbose = F, log.p = T)    
    res_log2 = relabelMCMC(res_log$permuted, maxit = 100, verbose = F, log.p = T)
    
    expect_equal(res_log2$iterations, 0L)
    expect_equal(res_log2$status, 0L)
    expect_equal(res_log2$permuted, res_log$permuted)
    expect_equal(res_log2$perms, matrix(rep(1:K, S), nr = S, byrow = TRUE))

})


test_that("relabelTRUE throws appropriate errors", {
    
    # params
    N = sample(20:100, 1)
    K = sample(2:4, 1)
    pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
    S = 100

    # generate array 
    test_ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    
    # sample true prob mat
    true = rdirichlet(N, pvec)
    
    # throw error if probs don't sum to one (due to dim mismatch)
    expect_error(relabelTRUE(test_ar, true, verbose = F, log.p = F))
    expect_error(relabelTRUE(aperm(test_ar, c(3,1,2)), t(true), F, F))
    
    # same test for log-probs
    expect_error(relabelTRUE(log(test_ar), log(true), F, T))
    expect_error(
      relabelTRUE(log(aperm(test_ar, c(3,1,2))), log(t(true)), F, T)
    )
    
    # randomly relabel some of the draws
    test_ar = aperm(test_ar, c(3,1,2))
    rel_draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel_draws) 
        test_ar[s, , ] = test_ar[s, , sample.int(K, K, F)]

    # check that it works for both verbose and silent
    expect_output(relabelTRUE(test_ar, true, T, F), regexp = ".")
    expect_silent(relabelTRUE(test_ar, true, F, F))
    expect_output(relabelTRUE(log(test_ar), log(true), T, T), regexp = ".")
    expect_silent(relabelTRUE(log(test_ar), log(true), F, T))
    
    # relabel
    res1 = relabelTRUE(test_ar, true, F, F)
    
    # relabel a second time
    test_ar2 = res1$permuted
    res2 = relabelTRUE(test_ar2, true, F, F)

    # first and second relabeling should be equal
    expect_equal(res1$permuted, res2$permuted)
    # but permutations should be different
    expect_false(isTRUE(all.equal(res1$perms, res2$perms)))
    # and perms for second permutation should be unchanged
    expect_equal(res2$perms, matrix(rep(1:K, S), nr = S, byrow = T))
    
    # relabel log
    res_log = relabelTRUE(log(test_ar), log(true), F, T) 
    
    # compare prob with log-prob results
    expect_false(isTRUE(all.equal(res1$permuted, res_log$permuted)))
    expect_equal(res1$permuted, exp(res_log$permuted))
    # compare perms
    expect_equal(res1$perms, res_log$perms)

})

test_that("relabling based on log-probs gives same results", {
    
    # how many tests to run
    n_test = 5
    
    for (tt in 1:n_test) {
        
        # params
        N = sample(20:100, 1)
        K = sample(2:4, 1)
        pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
        S = 100

        # generate array 
        test_ar = simplify2array(
                    lapply(1:S, function(w) rdirichlet(N, pvec))
                  )
        # reshuffle dimensions to S * N * K
        test_ar = aperm(test_ar, c(3, 1, 2))
        # randomly relable some of the draws
        rel_draws = sample.int(S, floor(S/4), replace = F)
        for (s in rel_draws) 
            test_ar[s, , ] = test_ar[s, , sample.int(K, K, F)]

        # results of this package
        res = relabelMCMC(test_ar, 100, verbose = F, log.p = F)
        res_log = relabelMCMC(log(test_ar), 100, verbose = F, log.p = T)
        
        # check dimensions
        expect_true(identical(dim(res_log$permuted), dim(test_ar)))
        # each label appears only once in each row
        expect_equal(
            sum(sapply(1:S, function(s) {
                length(res_log$perms[s,]) != length(unique(res_log$perms[s,]))
            })), 
            0
        )
        
        # compare results
        expect_equal(res_log$permuted, log(res$permuted))
        expect_true(sum(res$perms!=res_log$perms) == 0)
        
    }
  
    
})
