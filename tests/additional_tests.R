# Additional tests, not to be run on Cran or when testing the rest

library(rstan)
library(here)
library(testthat)
library(gridExtra)
library(relabelKL)

par.list = list(pi = c('pi', 'cols'),
                theta = c('theta', 'both'))


fit = readRDS(here("MMSBM.rds"))
re.fit = permute.stanfit(fit, par.list, log.form = F, maxit = 100, verbose = T)

# g1 = traceplot(fit, 'theta')
# g2 = traceplot(re.fit, 'theta')
# grid.arrange(grobs = list(g1, g2), nrow = 2)
# 
# 
# g1 = lapply(1:3, function(k) traceplot(fit, paste0('pi[1,',k,']')))
# g2 = lapply(1:3, function(k) traceplot(re.fit, paste0('pi[1,',k,']')))
# grid.arrange(grobs = c(g1, g2), nrow = 2)

re2.fit = permute.stanfit(re.fit, par.list, log.form = F, maxit = 100, verbose = T)

identical(re.fit, re2.fit)
