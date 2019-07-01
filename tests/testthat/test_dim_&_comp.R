context("Check Dimensions and Compare with label.swithching package")

has.ls = require(label.switching)

# function to draw from dirichlet distribution
rdirichlet = function(N, alpha) {
    
    l <- length(alpha)
    x <- matrix(rgamma(l * N, alpha), ncol = l, byrow = T)
    x/as.vector(x %*% rep(1, l))
}

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

        # randomly relable some of the draws
        rel.draws = sample.int(S, floor(S/4), replace = F)
        for (s in rel.draws) 
            test.ar[,,s] = test.ar[,sample.int(K, K, F) ,s]

        # results from the label.switching package
        if (has.ls) 
            res1 = stephens(aperm(test.ar, c(3,1,2)))
        # results of this package
        res2 = relabelMCMC(test.ar, 100, verbose = F)
        
        # check dimensions
        expect_true(identical(dim(res2$relabeled), dim(test.ar)))
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

    test.ar = simplify2array(
                    lapply(1:S, function(w) rdirichlet(N, pvec))
                  )

    # randomly relable some of the draws
    rel.draws = sample.int(S, floor(S/4), replace = F)
    for (s in rel.draws) 
        test.ar[,,s] = test.ar[,sample.int(K, K, F) ,s]

    res = relabelMCMC(test.ar, verbose = F)    
    
    expect_equal(permuteMCMC(test.ar, res$perms), res$relabeled)
    
    samp.ind = sample.int(dim(test.ar)[1], 1L)
    expect_equal(
        permuteMCMC(t(test.ar[samp.ind,,]), res$perms),
        t(res$relabeled[samp.ind,,])
    )

})