# #' Update theta estimate
# #' @import spam
# 
# VB_update_theta <- function(state) {
#   ## helper function:
#   lambda <- function(x)  -tanh(x/2)/(4*x)
#   #'
#   state$m_old <- m <- state$m
#   S <- state$S
#   #' Constants
#   p <- state$z
#   Phip <- state$Phi %*% p
#   Phi2 <- state$Phi^2
#   pPhi2ppb <- p*Phi2%*%(p-p^2)
#   v <- (Phip)^2 + Phi2%*%(p*(1-p))
#   X <- cbind.spam(1-p, p, p*Phip)
#   oo <- state$offset^2
#   #' Initial:
#   xi <- state$xi
#   #'
#   loop <- TRUE
#   while(loop) {
#     xi_old <- xi
#     la <- lambda(xi)
#     L <- diag.spam(x=la)
#     #' Update S:
#     A <- sum(X[,1] * la)
#     B <- sum(X[,2] * la)
#     C <- sum(X[,3] * la)
#     D <- (t(X[,3])%*%L%*%Phip)[1] + sum(pPhi2ppb * la)
#     H <- matrix(c(A,0,0,0,B,C,0,C,D), ncol=3)
#     S <- solve( state$S0i - 2*H )
#     #' update m:
#     m <- S%*%t(t(state$y - 0.5 + 2 * L %*% state$offset)%*% X + state$Sm0)
#     #'
#     #' Update xi:
#     M <- S + m%*%t(m)
#     xi2 <- M[1,1]*(1-p) + 
#            M[2,2]*p +  
#            M[3,3]*p*v +
#          2*M[2,3]*p*Phip + 2 * diag( X%*%m%*%t(state$offset) ) + oo
#     xi <- sqrt(xi2)[,1]
#     #' check convergence
#     loop <- max( abs(xi - xi_old) ) > state$eps
#   }
#   
#   state$xi <- xi
#   state$m <- m
#   state$S <- S
#   state$m_hist <- cbind(state$m_hist, m[,1])
#   state$L <- L
#   state$la <- la
#   state
# }
# 
# 
# 
# #' Update z, more correct version
# #' @import spam
# VB_update_z <- function(state) {
#   #' Helper:
#   sig <- function(x) 1/(1+exp(-x))
#   #' updateable z
#   ind <- which(state$y == 1)
#   ord <- ind
#   #'
#   #' constants
#   L <- state$L
#   la <- diag(L)
#   m <- state$m[,1]
#   M <- state$S + m%*%t(m)
#   b <- state$y - 0.5  + 2 * L%*%state$offset
#   Phi <- state$Phi
#   Phi2 <- Phi^2
#   #'
#   state$z_old <- p <- state$z
#   #' update:
#   p0 <- p1 <- p
#   for(i in ord) {
#     p1[i] <- 1
#     p0[i] <- 0
#     #X1 <- cbind.spam(1-p1, p1, p1*Phi%*%p1)
#     #X0 <- cbind.spam(1-p0, p0, p0*Phi%*%p0)
#     #X <- cbind.spam(p0-p1, p1-p0, p1*state$Phi%*%p1 - p0*Phi%*%p0)
#     #F1 <- t(b)%*%(X1-X0)
#     #F1 <- t(b)%*%X
#     #V1 <- sum(m*F1)
#     #V1 <- m[1]*t(b)%*%(p0-p1)+m[2]*t(b)%*%(p1-p0) +
#     #      m[3]*t(b)%*%(p1*state$Phi%*%p1 - p0*Phi%*%p0)
#     spp <- sum(Phi[i,]*p)
#     ppj <- Phi[,i]*p
#     V1 <- -m[1]*b[i]+m[2]*b[i] + m[3]*( sum(b*ppj)+b[i]*spp ) 
#     #
#     A <- -L[i,i]
#     B <- -A
#     #C <- L[i,i]*Phi[i,]%*%p + t(Phi[,i])%*%L%*%p
#     C <- B*spp + sum(la*ppj)
#     #D1 <- t(p1)%*%L%*%((Phi%*%p1)^2+Phi2%*%(p1*(1-p1)))
#     D1 <- sum( p1*la*((Phi%*%p1)^2+Phi2%*%(p1*(1-p1))) )
#     #D0 <- t(p0)%*%L%*%((Phi%*%p0)^2+Phi2%*%(p0*(1-p0)))
#     D0 <- sum( p0*la*((Phi%*%p0)^2+Phi2%*%(p0*(1-p0))) )
#     D <- D1-D0
#     V2 <- M[1,1]*A + M[2,2]*B + M[3,3]*D + 2*M[2,3]*C
#     u <- V1 + V2 + state$p00
#     p0[i]<-p1[i]<-p[i] <- sig(u[1])
#   }
#   state$z[ind] <- p[ind]
#   state$z_hist <- rbind(state$z_hist, p[ind])
#   state  
# }
# 
# 
# 
# #' Update z estimate, rough mean estimate
# VB_update_z_rough <- function(state) {
#   #' Helper:
#   sig <- function(x) 1/(1+exp(-x))
#   #' updateable z
#   ind <- which(state$y == 1)
#   ord <- ind
#   #'
#   #' constants
#   L <- state$L
#   m <- state$m[,1]
#   M <- state$S + m%*%t(m)
#   b <- state$y - 0.5  + 2 * L%*%state$offset
#   #'
#   state$z_old <- p <- state$z
#   #' update:
#   for(i in ord) {
#     p0 <- p1 <- p
#     p1[i] <- 1
#     p0[i] <- 0
#     X1 <- cbind.spam(1-p1, p1, p1*state$Phi%*%p1)
#     X0 <- cbind.spam(1-p0, p0, p0*state$Phi%*%p0)
#     F1 <- t(b)%*%(X1-X0)
#     V1 <- sum(m*F1)
#     H <- t(X1)%*%L%*%X1 - t(X0)%*%L%*%X0
#     V2 <- sum( diag( M%*%H  ) )
#     u <- V1 + V2 + state$p00
#     p[i] <- sig(u[1])
#   }
#   state$z[ind] <- p[ind]
#   state$z_hist <- rbind(state$z_hist, p[ind])
#   state  
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #' Check VB convergence
# VB_check_convergence <- function(state){
#   dig <- 4
#   mdelta <- round(abs( state$m - state$m_old  ), dig)
#   zdelta <- round(max( abs(state$z - state$z_old) ), dig)[1]
#   state$loop <- max(mdelta, zdelta) > state$eps
#   state$diag <- paste0("[m ", paste0(mdelta[,1], collapse=","), "] [z ", zdelta,"]")
#   state
# }
# 
# 
# 
# 
# 
# 
