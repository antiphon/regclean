#' test alt version

library(devtools)
load_all(".")


tt <- system.time( ff <- MCMC_ApproximationAlphaClassify_alt(obs, R, mcmcPar = list(reps=1e4, repin=1e2)) )

print(tt)


tt <- system.time( ff <- VBclassify(obs, R, verb=T, prior=list(S=diag(c(100,100,10)))))

print(tt)


