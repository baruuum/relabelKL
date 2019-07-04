context("Permuting and monitoring stanfit objects")


test_that("extract and back-transformation works", {
    
    skip_on_cran()
    skip_on_travis()
    
    data("mmsbm")
    
    stan.pi = extract_n_combine(mmsbm, "pi")

    expect_true(
        isTRUE(
            all.equal(
                to_stan_array(stan.pi),
                rstan:::extract(mmsbm, par = 'pi', permute = F, inc_warmup = F)
            )
        )
    )
    
    stan.theta = extract_n_combine(mmsbm, "theta")
    arr.theta = to_stan_array(stan.theta)
    
    expect_true(
        isTRUE(
            all.equal(
                arr.theta,
                rstan:::extract(mmsbm, par = 'theta', permute = F, inc_warmup = F)
            )
        )
    )
    
    skip("skip producing plots")
    array_traceplot(arr.theta, 'theta')
    array_traceplot(arr.theta, 'theta[1,1]')
    
    rel = relabelMCMC(stan.pi, 50, TRUE)
    re.theta = permuteMCMC(stan.theta, rel$perms, "both")
    array_traceplot(to_stan_array(re.theta), "theta")
    
})

