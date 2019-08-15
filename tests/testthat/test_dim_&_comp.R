context("Check Dimensions and Compare with label.swithching package")

has.ls = requireNamespace("label.switching", quietly = T)

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
    test.ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    # reshuffle dimensions to match arg
    test.ar = aperm(test.ar, c(3, 1, 2))
    
    expect_error(relabelMCMC(test.ar, 100, F, T))
    expect_error(relabelMCMC(log(test.ar), 100, F, F))
  
})

test_that("results coincide with stephens function from the label.switching package", {
    
    # how many tests to run
    n.test = 5
    
    for (tt in 1:n.test) {
        
        # params
        N = sample(20:100, 1)
        K = sample(2:4, 1)
        pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
        S = 100

        # generate array 
        test.ar = simplify2array(
                    lapply(1:S, function(w) rdirichlet(N, pvec))
                  )
        # reshuffle dimensions to S * N * K
        test.ar = aperm(test.ar, c(3, 1, 2))
        # randomly relable some of the draws
        rel.draws = sample.int(S, floor(S/4), replace = F)
        for (s in rel.draws) 
            test.ar[s, , ] = test.ar[s, , sample.int(K, K, F)]

        # results from the label.switching package
        if (has.ls) 
            res1 = label.switching::stephens(test.ar)
        
        # results of package
        res2 = relabelMCMC(test.ar, maxit = 100, verbose = F, log.p = F)
        
        # check dimensions
        expect_true(identical(dim(res2$permuted), dim(test.ar)))
        # each label appears only once in each row
        expect_equal(
            sum(sapply(1:S, function(s) {
                length(res2$perms[s,]) != length(unique(res2$perms[s,]))
            })), 
            0
        )
        
        # compare with label.switching::stephens
        if (has.ls)
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
    test.ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    # reshuffle dimensions to S * N * K
    test.ar = aperm(test.ar, c(3, 1, 2))
    # randomly relabel some of the draws
    rel.draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel.draws) 
        test.ar[s, , ] = test.ar[s, , sample.int(K, K, F)]
    
    expect_error(
        relabelMCMC(aperm(test.ar, c(2,3,1)), verbose = F, log.p = F)
    )

    res = relabelMCMC(test.ar, verbose = F, log.p = F)    
    
    expect_equal(permuteMCMC(test.ar, res$perms, "cols"), res$permuted)
    expect_error(permuteMCMC(test.ar, res$perms, "both"))
    
    samp.ind = sample.int(dim(test.ar)[2], 1L)
    expect_equal(
        permuteMCMC(test.ar[, samp.ind, ], res$perms, "cols"),
        res$permuted[, samp.ind,]
    )
    expect_error(permuteMCMC(test.ar[, samp.ind, ], res$perms, "rows"))
    expect_error(permuteMCMC(test.ar[, samp.ind, ], res$perms, "both"))
    
})

test_that("2nd relabeling results in identity mapping", {
    
    # params
    N = sample(20:100, 1)
    K = sample(2:4, 1)
    pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
    S = 100

    # generate array 
    test.ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    # reshuffle dimensions to S * N * K
    test.ar = aperm(test.ar, c(3, 1, 2))
    # randomly relabel some of the draws
    rel.draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel.draws) 
        test.ar[s, , ] = test.ar[s, , sample.int(K, K, F)]

    res = relabelMCMC(test.ar, maxit = 100, verbose = F, log.p = F)    
    res2 = relabelMCMC(res$permuted, maxit = 100, verbose = F, log.p = F)
    
    expect_equal(res2$iterations, 0L)
    expect_equal(res2$status, 0L)
    expect_equal(res2$permuted, res$permuted)
    expect_equal(res2$perms, matrix(rep(1:K, S), nr = S, byrow = TRUE))
    
    # same test for log-scale
    res.log = relabelMCMC(log(test.ar), maxit = 100, verbose = F, log.p = T)    
    res.log2 = relabelMCMC(res.log$permuted, maxit = 100, verbose = F, log.p = T)
    
    expect_equal(res.log2$iterations, 0L)
    expect_equal(res.log2$status, 0L)
    expect_equal(res.log2$permuted, res.log$permuted)
    expect_equal(res.log2$perms, matrix(rep(1:K, S), nr = S, byrow = TRUE))

})

test_that("relabelTRUE throws appropriate errors", {
    
    # params
    N = sample(20:100, 1)
    K = sample(2:4, 1)
    pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
    S = 100

    # generate array 
    test.ar = simplify2array(
                lapply(1:S, function(w) rdirichlet(N, pvec))
              )
    
    # sample true prob mat
    true = rdirichlet(N, pvec)
    
    # throw error if probs don't sum to one
    expect_error(relabelTRUE(test.ar, true, verbose = F, log.p = F))
    expect_error(relabelTRUE(aperm(test.ar, c(3,1,2)), t(true), F, F))
    
    # same test for log-probs
    expect_error(relabelTRUE(log(test.ar), log(true), F, T))
    expect_error(
      relabelTRUE(log(aperm(test.ar, c(3,1,2))), log(t(true)), F, T)
    )
    
    # randomly relabel some of the draws
    test.ar = aperm(test.ar, c(3,1,2))
    rel.draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel.draws) 
        test.ar[s, , ] = test.ar[s, , sample.int(K, K, F)]

    
    # relabel
    res1 = relabelTRUE(test.ar, true, F, F)
    
    # relabel a second time
    test.ar2 = res1$permuted
    res2 = relabelTRUE(test.ar2, true, F, F)

    # first and second relabeling should be equal
    expect_equal(res1$permuted, res2$permuted)
    # but permutations should be different
    expect_false(isTRUE(all.equal(res1$perms, res2$perms)))
    # and perms for second permutation should be unchanged
    expect_equal(res2$perms, matrix(rep(1:K, S), nr = S, byrow = T))
    
    # relabel log
    res.log = relabelTRUE(log(test.ar), log(true), F, T) 
    
    # compare prob with log-prob results
    expect_false(isTRUE(all.equal(res1$permuted, res.log$permuted)))
    expect_equal(res1$permuted, exp(res.log$permuted))
    # compare perms
    expect_equal(res1$perms, res.log$perms)

})

test_that("relabling based on log-probs gives same results", {
    
    # how many tests to run
    n.test = 5
    
    for (tt in 1:n.test) {
        
        # params
        N = sample(20:100, 1)
        K = sample(2:4, 1)
        pvec = runif(K, 0.1, 1.0) * runif(1, 0.1, 10)
        S = 100

        # generate array 
        test.ar = simplify2array(
                    lapply(1:S, function(w) rdirichlet(N, pvec))
                  )
        # reshuffle dimensions to S * N * K
        test.ar = aperm(test.ar, c(3, 1, 2))
        # randomly relable some of the draws
        rel.draws = sample.int(S, floor(S/4), replace = F)
        for (s in rel.draws) 
            test.ar[s, , ] = test.ar[s, , sample.int(K, K, F)]

        # results of this package
        res = relabelMCMC(test.ar, 100, verbose = F, log.p = F)
        res.log = relabelMCMC(log(test.ar), 100, verbose = F, log.p = T)
        
        # check dimensions
        expect_true(identical(dim(res.log$permuted), dim(test.ar)))
        # each label appears only once in each row
        expect_equal(
            sum(sapply(1:S, function(s) {
                length(res.log$perms[s,]) != length(unique(res.log$perms[s,]))
            })), 
            0
        )
        
        # compare results
        expect_equal(res.log$permuted, log(res$permuted))
        expect_true(sum(res$perms!=res.log$perms) == 0)
        
    }
  
    
})
