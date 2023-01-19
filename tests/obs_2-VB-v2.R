#' VB, test 1
library(rstrauss)
library(spatstat)
library(devtools)
load_all(".")


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
lambda <- intensity(x0)
obs <- superimpose(x1,x0)
truth <- rep(c(1,0), c(n1,n0))

#' subset, as returned by the MCMC
sub <- inside.owin(coords(obs),w = W0)
obss <- obs[sub]
truths <- truth[sub]
obss_x <- list(x=cbind(obss$x, obss$y), bbox=cbind(0:1, 0:1))

#'
#' Then classify with VB
#' 
fit <- VBclassify(R=R, x=obs, v2=FALSE)
fit2 <- VBclassify(R=R, x=obs, v2=TRUE)

#'
#'
comp <- data.frame(true=c(lambda,beta,gamma), est=fit$coef, est2=fit2$coef)
print(round(comp, 3))
print(c(fit$took ,fit2$took))
#' Diagnose:
par(mfrow=c(2,2))
plot(fit$prob, ylim=c(0,1), col=truth+2)
plot(fit$z_hist[,1], ylim=c(0,1), type="l")
for(i in 1:ncol(fit$z_hist))lines(fit$z_hist[,i], col=i)

plot(fit2$prob, ylim=c(0,1), col=truth+2)
plot(fit2$z_hist[,1], ylim=c(0,1), type="l")
for(i in 1:ncol(fit2$z_hist))lines(fit2$z_hist[,i], col=i)


