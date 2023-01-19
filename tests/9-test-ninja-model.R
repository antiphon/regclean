#' Test the parameter free model, working flag 'ninja'

library(devtools)
load_all(".")

load("test_pat1.rda")

x <- VB_check_pattern(obs)

z <- VB_check_prior( x = x, ninja = TRUE)

f <- VBclassify(x, R =R , verb =T)

f2 <- VBclassify(x, R=R, ninja =T, verb =T )


par(mfrow=c(2,2))
plot(f$prob, ylim=c(0,1), col=truth+2)
plot(f)

plot(f2$prob, ylim=c(0,1), col=truth+2)
plot(f2)
