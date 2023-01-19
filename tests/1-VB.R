#' VB, test  v2, different update for z_i
#library(spam)
library(rstrauss)
library(spatstat)
library(devtools)
load_all(".")

set.seed(1)
#' data
R <- 0.04
W <- as.owin(c(-R,1+R,-R,1+R))
W0 <- square()
gamma <- 0.01
beta <- 200

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
#'
#' Then classify with VB
#' 
fit <- VBclassify(R=R, x=obs, verb=T)
#'
#'
comp <- data.frame(true=c(lambda,beta,gamma), est=fit$coef)
print(round(comp, 3))
#' Diagnose:
par(mfrow=c(2,2))
plot(fit$prob, ylim=c(0,1), col=truth+2)
plot(fit$z_hist[,1], ylim=c(0,1), type="l")
for(i in 1:ncol(fit$z_hist))lines(fit$z_hist[,i], col=i)


summary(fit)
