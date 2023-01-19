#' Test VB for large pattern

library(spam)
library(rstrauss)
library(devtools)
load_all(".")

set.seed(2)
#' data
gamma <- .1
beta <- 1500
R <- strauss_theory(beta, 2)*0.8
noiselevel <- 0.2



#' true
x1 <- rstrauss(beta, gamma, R, bbox=cbind(0:1, 0:1), perfect=TRUE,  iter = 2e5, verb=!T)
cat("sim'd\n")
n1 <- nrow(x1$x)
print(n1)
#' noise:
n0 <- round(  n1*noiselevel/(1-noiselevel)  )
x0 <- matrix(runif(n0*2), ncol=2)
lambda <- n0
obs <- list(x=rbind(x1$x, x0), bbox=x1$bbox)
truth <- rep(c(1,0), c(n1,n0))
#'
#' Then classify with VB
#' 
#Rprof()
set.seed(1)
fit <- VBclassify(R=R, x=obs, verb=T)
#s<-summaryRprof()
#print(s)
print(fit$took)

print(fit$coef)
ss <- summary(fit)
plot(fit$prob, ylim=0:1, col=1+truth, pch=19)
plot(fit)
