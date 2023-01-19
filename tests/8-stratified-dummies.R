#' stratified sampling for dummies
#' 

library(devtools)
load_all(".")

load("test_pat1.rda")

x <- VB_check_pattern(obs)

prior <- VB_check_prior(NULL,x)


z <- VB_initial_state(x, R, prior, quads = 20, eps = 1e-3, keep.Q = T, rho=90)
