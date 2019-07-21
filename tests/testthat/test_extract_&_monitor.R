context("Permuting and monitoring stanfit objects")

test_that("extract and back-transformation works", {
    
    data("mmsbm")
    
    # check whether extract_n_combine -> to_stan_array returns original array
    stan.pi = extract_n_combine(mmsbm, "pi")
    ext.pi = rstan:::extract(mmsbm, par ="pi", permute = F)
    expect_true(identical(ext.pi, to_stan_array(stan.pi)))
    
    stan.theta = extract_n_combine(mmsbm, "theta")
    ext.theta = rstan:::extract(mmsbm, par = "theta", permute = F)
    expect_true(identical(ext.theta, to_stan_array(stan.theta)))
    
    skip_on_cran()
    skip_on_travis()
    
    arr.theta = to_stan_array(stan.theta)
    
    skip("skip producing plots")
    array_traceplot(arr.theta, 'theta')
    array_traceplot(arr.theta, 'theta[1,1]')
    
    rel = relabelMCMC(stan.pi, 50, TRUE, FALSE)
    re.theta = permuteMCMC(stan.theta, rel$perms, "both")
    array_traceplot(to_stan_array(re.theta), "theta")
    
})

