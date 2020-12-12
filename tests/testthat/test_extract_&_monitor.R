context("Permuting and monitoring stanfit objects")

# skip on 32bit windows (have to figure out what's going wrong here)
# if (.Platform$OS.type == "windows" && .Platform$r_arch == "i386")
#     skip_on_appveyor()

test_that("extract and back-transformation works", {
    
    data("mmsbm")
    
    # check whether non-stanfit object results in error
    a = unclass(mmsbm)
    expect_error(extract_n_combine(a, "pi"))
    
    # check whether non-character second argument throws error
    expect_error(extract_n_combine(mmsbm, pi))
    
    # check whether extract_n_combine -> to_stan_array returns original array
    expect_true(
        identical(
            rstan:::extract(mmsbm, par ="pi", permute = F), 
            to_stan_array(extract_n_combine(mmsbm, "pi"))
        )
    )
    
    expect_true(
        identical(
            rstan:::extract(mmsbm, par = "theta", permute = F), 
            to_stan_array(extract_n_combine(mmsbm, "theta"))
        )
    )
    
    # check whether plots are produced by array_traceplot
    
    arr_theta = to_stan_array(extract_n_combine(mmsbm, "theta"))
    
    a = array_traceplot(arr_theta, "theta")
    expect_equal(class(a)[2], "ggplot")
    a = array_traceplot(arr_theta, "theta[1,1]")
    expect_equal(class(a)[2], "ggplot")
    a = array_traceplot(arr_theta, c("theta[1,1]", "theta[2,2]"))
    expect_equal(class(a)[2], "ggplot")
    
    expect_error(array_traceplot(arr_theta, c("beta")))
    expect_error(array_traceplot(arr_theta, c("theta[4,2]")))
    
    # check whether relabeling routine works
    stan_pi = extract_n_combine(mmsbm, "pi")
    stan_theta = extract_n_combine(mmsbm, "theta")
    rel = relabelMCMC(stan_pi, log_p = F, renormalize = F, 50, 0, FALSE)
    re_theta = permuteMCMC(stan_theta, rel$perms, "both")
    a = array_traceplot(to_stan_array(re_theta), "theta")
    expect_equal(class(a)[2], "ggplot")
    
    # check whether rstan::monitor function can be applied
    mres = rstan::monitor(to_stan_array(stan_theta))
    expect_equal(attr(stan_theta, "org.att")$dimnames$parameters,
                 rownames(mres))
    
})

