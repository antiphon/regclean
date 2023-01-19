#' optimize vb, v1
#' 

library(devtools)
load_all(".")

load("test_pat1.rda")

 x<- cbind(obs$x, obs$y)

set.seed(1)
t1 <- system.time( e1 <- VBclassify(x, R=R, verb=TRUE, v2=TRUE) )
set.seed(1)
t0 <- system.time( e0 <- VBclassify(x, R=R, verb=TRUE, v2=FALSE) )

print(rbind(t0,t1))

print(cbind(e0$coef, e1$coef))
