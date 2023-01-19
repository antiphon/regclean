#' Noise removal from regular point patterns
#' 
#' Base process assumed to be the Strauss model, noise process assumed Poisson, processes independent.
#' 
#' @param x Point pattern, as list with $x coordinate matrix $bbox bound box ranges (col-wise) (or a ppp)
#' @param R Range for Strauss model.
#' @param prior Priors, see Details.
#' @param eps Convergence criterion.
#' @param rho Intensity of dummies, default = 2*n/Volume. Two dummy sets will be generated.
#' @param verb Verbose output?
#' @param quads How many subwindows to use for dummy sampling stratification? Default is 1, no stratification.
#' @param ninja Instead of Strauss potential gamma^1(d < R), use a diminishing potential without extra parameter: (d/R)^1(d < R)
#' @details
#' Only rectangle windows supported at the moment.
#' 
#' The priors: It should be a list.
#'
#' @references 
#' [1] A. Baddeley, J. Coeurjolly, E. Rubak, and R. Waagepetersen, "Logistic regression for spatial Gibbs point processes," Biometrika, 2013.
#' 
#' @export

VBclassify <- function(x, R, prior, eps=0.01, rho, verb=FALSE, quads = 1, ninja = FALSE) {
  T0 <- Sys.time()
  #
  # Verbose?
  cat2 <- if(verb) cat else function(...) NULL
  # Check pattern
  x <- VB_check_pattern(x)
  # Check priors
  prior <- VB_check_prior(prior, x, ninja = ninja)
  #
  # Create the state-object
  state <- VB_initial_state(x, R, prior, rho, eps, quads, ninja = ninja)
  #
  #
  #
  # Main loop:
  iter <- 0
  loop <- TRUE
  while(loop) {
    # Update theta
    state <- VB_update_theta(state)
    # Update z
    state <- VB_update_z(state)
    # check convergence
    state <- VB_check_convergence(state)
    #
    loop <- state$loop
    iter <- iter + 1
    cat2(state$diag, "    \n")
  }
  cat2("\n")
  #
  # Format result object:
  # 
  post <- list(m=state$m, S=state$S)
  prob <- state$z[1:nrow(x$x)]
  pred <- 1 * (prob>0.5)
  z_hist <- state$z_hist
  # mean estimates in typical parametrization (non exp-family)
  coef <- exp(post$m[,1] + diag(post$S)/2)
  names(coef) <- c("lambda", "beta", "gamma")
  #
  # compile:
  out <- list(coef=coef,
              pred=pred,
              prob=prob,
              posterior=post,
              z_hist = z_hist,
              call = sys.call(),
              prior=prior,
              R=R,
              eps=eps,
              state=state,
              iter = iter, xi_iter = state$xi_iter,
              data_facts = list(n=nrow(x$x), bbox=x$bbox),
              took=Sys.time()-T0)
  class(out) <- c("VBclassification", is(out))
  # Done!
  out
}
