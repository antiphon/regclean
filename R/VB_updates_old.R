# #' Update theta estimate
# #' @import Matrix
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
#   X <- cBind(1-p, p, p*Phip)
#   oo <- state$offset^2
#   #' Initial:
#   xi <- state$xi
#   #'
#   loop <- TRUE
#   while(loop) {
#     xi_old <- xi
#     la <- lambda(xi)
#     L <- Diagonal(x=la)
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
#   state
# }
# 
# 
# 
# #' Update z estimate
# 
# VB_update_z <- function(state) {
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
#     X1 <- cBind(a=1-p1, b=p1, c=p1*state$Phi%*%p1)
#     X0 <- cBind(a=1-p0, b=p0, c=p0*state$Phi%*%p0)
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
# #' Check VB convergence
# VB_check_convergence <- function(state){
#   dig <- 4
#   mdelta <- round(abs( state$m - state$m_old  ), dig)
#   zdelta <- round(max( abs(state$z - state$z_old) ), dig)
#   state$loop <- max(mdelta, zdelta) > state$eps
#   state$diag <- paste0("[m ", paste0(mdelta, collapse=","), "] [z ", zdelta,"]")
#   state
# }
# 
# 
# 
# 
# 
# 
