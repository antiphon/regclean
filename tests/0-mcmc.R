#' MCMC, test 1
library(rstrauss)
library(regclean)
library(devtools)
load_all(".")

lambert_W0 <- LambertW # something wrong, need to define this.


#' data
R <- 0.08
W <- as.owin(c(-R,1+R,-R,1+R))
W0 <- square()
gamma <- 0.01
beta <- 100

noiselevel <- 0.2

#' true
x1 <- rStrauss(beta = beta, gamma = gamma, R = R, W=W)
n1 <- x1$n
#' noise:
n0 <- round(  n1*noiselevel/(1-noiselevel)  )
x0 <- runifpoint(n0, win = W)

obs <- superimpose(x1,x0)
truth <- rep(c(1,0), c(n1,n0))
#' subset, as returned by the MCMC
sub <- inside.owin(coords(obs),w = W0)
obss <- obs[sub]
truths <- truth[sub]
obss_x <- list(x=cbind(obss$x, obss$y), bbox=cbind(0:1, 0:1))


#' Then try MCMC
#' 
fit <- classificationMCMCApproximationAlpha(R, obs, 1e4, 1e3)
#' check also VB
fitvb <- VBclassify(obs, R)
#
par(mfrow=c(3,1))
p <- function(x) plot(x, legend = F, cols=1:2+2, cex=1.2)
p(setmarks(obs, factor(truth) )); plot(W0, add=T)
p(setmarks(obss, factor(1*(fit$p>0.5), levels=0:1) ))
p(setmarks(obs, factor(fitvb$pred ))); plot(W0, add=T)

#' compare
comp <- data.frame(true=c(lambda=n0/area.owin(W),beta=beta, gamma=gamma), mcmc=fit$theta[-4], vb=fitvb$coef)
print(comp)


