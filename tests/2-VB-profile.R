#' profile v2

library(Matrix)
library(spam)
library(rstrauss)
library(spatstat)
library(devtools)
load_all(".")

set.seed(2)
#' data
R <- 0.08
W <- as.owin(c(-R,1+R,-R,1+R))
W0 <- square()
gamma <- 0.01
beta <- 50

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
#Rprof()
set.seed(1)
fit <- VBclassify(R=R, x=obs, verb=!T)
#s<-summaryRprof()
print(fit$took)

print(fit$coef)
