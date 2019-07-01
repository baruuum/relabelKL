context("Check Dimensions and Compare with label.swithching package")

has.ls = require(label.switching)

# function to draw from dirichlet distribution
rdirichlet = function(N, alpha) {
    
    l <- length(alpha)
    x <- matrix(rgamma(l * N, alpha), nc = l, byrow = T)
    x/as.vector(x %*% rep(1, l))
}

test_that("results coincide with stephens function from the label.switching package", {
    
    # params
    N = 50
    pvec = c(2,3,1) * 0.5
    K = length(pvec)
    S = 100
    
    # how many tests to run
    n.test = 5
    
    for (tt in 1:n.test) {

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
        res2 = relabel(test.ar, verbose = F)
        
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
