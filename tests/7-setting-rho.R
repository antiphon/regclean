# test rho setup
#' VB, test 1
library(devtools)
load_all(".")

load("test_pat1.rda")

#' Classify with VB: how does impalance in rho affect?
#' 
fit <- VBclassify(R=R, x=obs, rho=200, verb=TRUE)
fit2 <- VBclassify(R=R, x=obs, rho = c(100, 200), verb=TRUE)

#'
#'
comp <- data.frame( est=fit$coef, est2=fit2$coef)
print( round(comp, digits=3))
print(c(fit$took ,fit2$took))
#' Diagnose:
par(mfrow=c(2,2))
plot(fit$prob, ylim=c(0,1), col=truth+2)
plot(fit)

plot(fit2$prob, ylim=c(0,1), col=truth+2)
plot(fit2)


