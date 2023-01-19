#' Profile mcmc v1
library(rstrauss)
library(spatstat)

library(devtools)
load_all(".")
library(coda)

lambert_W0 <- LambertW

load("test_pat1.rda")

Rprof()
fit <- MCMC_ApproximationAlphaClassify(obs, R, mcmcPar = list(reps=1e3, repin=1e2))
s <- summaryRprof()
