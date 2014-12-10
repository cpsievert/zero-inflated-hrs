metropolis <- function(beta, alpha, delta, lambda, jump_sd, param="beta") {
  t <- exp(rep(beta, ni) + rep(alpha, ni)*dat$home)  #save repeated calc for theta (p is for potential)
  theta <- t/(1 + t)
  p <- exp(rep(delta, ni) + rep(lambda, ni)*dat$home)  #save repeated calc for lambda
  psi <- p/(1 + p)
  loglike <- ifelse(zero, log((1-theta) + theta*(1-psi)^dat$n), 
                    dat$y*log(psi) + (dat$n - dat$y)*log(1-psi) + log(theta))
  #compute likelihood at the player level
  denom <- diff(c(0, cumsum(loglike)[cum.ni]))
  
  if (param == "beta") {
    prior.old <- dnorm(beta, 0, 5, log=TRUE) #controls prior belief?
    beta_star <- rnorm(n.players, mean=beta, sd=jump_sd)
    if(any(!is.finite(beta_star))) browser()
    prior.star <- dnorm(beta_star, 0, 5, log=TRUE)
    t <- exp(rep(beta_star, ni) + rep(alpha, ni)*dat$home)  
    theta <- t/(1 + t)
  } 
  if (param == "alpha") {
    prior.old <- dnorm(alpha, 0, 5, log=TRUE)
    alpha_star <- rnorm(n.players, mean=alpha, sd=jump_sd)
    if(any(!is.finite(alpha_star))) browser()
    prior.star <- dnorm(alpha_star, 0, 5, log=TRUE)
    t <- exp(rep(beta, ni) + rep(alpha_star, ni)*dat$home)
    theta <- t/(1 + t)
  } 
  if (param == "delta") {
    prior.old <- dnorm(delta, 0, 5, log=TRUE) #sigma_delta makes prior diffuse
    delta_star <- rnorm(n.players, mean=delta, sd=jump_sd)
    if(any(!is.finite(delta_star))) browser()
    prior.star <- dnorm(delta_star, 0, 5, log=TRUE)
    p <- exp(rep(delta_star, ni) + rep(lambda, ni)*dat$home)  
    psi <- p/(1 + p)
  } 
  if (param == "lambda") {
    prior.old <- dnorm(lambda, 0, 5, log=TRUE)
    lambda_star <- rnorm(n.players, mean=lambda, sd=jump_sd)
    if(any(!is.finite(lambda_star))) browser()
    prior.star <- dnorm(lambda_star, 0, 5, log=TRUE)
    p <- exp(rep(delta, ni) + rep(lambda_star, ni)*dat$home)
    psi <- p/(1 + p)
  } 
  #density in numerator of acceptance ratio
  num <- diff(c(0, cumsum(ifelse(zero, log((1-theta) + theta*(1-psi)^dat$n), 
                                 dat$y*log(psi) + (dat$n - dat$y)*log(1-psi) + log(theta)))[cum.ni]))
  ratio <- exp(num + prior.star - denom - prior.old)
  #lens <- c(length(num), length(denom), length(prior.star), length(prior.old))
  #cat(lens, "\n")
  #if (any(!is.finite(ratio))) browser()
  # impose acceptance probability
  if (param == "beta") {
    beta <- ifelse(runif(n.players) < ratio, beta_star, beta)
    return(list(beta=beta, ratio=ratio))
  } 
  if (param == "delta") {
    delta <- ifelse(runif(n.players) < ratio, delta_star, delta)
    return(list(delta=delta, ratio=ratio))
  } 
  if (param == "alpha") {
    alpha <- ifelse(runif(n.players) < ratio, alpha_star, alpha)
    return(list(alpha=alpha, ratio=ratio))
  } 
  if (param == "lambda") {
    lambda <- ifelse(runif(n.players) < ratio, lambda_star, lambda)
    return(list(lambda=lambda, ratio=ratio))
  } 
}

#MCMC workhorse function

#dat should be a data frame with:
# (1) 'url' - unique identifier for each game
# (2) 'batter_name' - unique identifier of players
# (3) 'home' - indicator for home vs. away game)
# (4) 'y' - number of occurences for outcome of interest
# (5) 'n' - number of atbats

mcmc <- function(dat, alpha.init = 0, lambda.init = 0, beta.init = 0, delta.init = 0,
                 n.reps = 3000, adapt = 1000, tune = TRUE){
  require(MASS)
  delta_keep <- beta_keep <- lambda_keep <- alpha_keep <- matrix(numeric(n.reps*n.players), 
                                                                 nrow=n.reps, ncol=n.players)
  #impose starting values
  alpha_keep[1,] <- alpha <- if (length(alpha.init)==1) rep(alpha.init, n.players) else alpha.init
  lambda_keep[1,] <- lambda <- if (length(lambda.init)==1) rep(lambda.init, n.players) else lambda.init
  beta_keep[1,] <- beta <- if (length(beta.init)==1) rep(beta.init, n.players) else beta.init
  delta_keep[1,] <- delta <- if (length(delta.init)==1) rep(delta.init, n.players) else delta.init
  
  #calculate log-likelihood for the starting value
  t <- exp(rep(beta, ni) + rep(alpha, ni)*dat$home)  #save repeated calc for theta
  theta <- t/(1 + t)
  p <- exp(rep(delta, ni) + rep(lambda, ni)*dat$home)  #save repeated calc for lambda
  psi <- p/(1 + p)
  #should I worry about numerical stability when dat$y==0???
  loglik <- sum(ifelse(zero, log((1-theta) + theta*(1-psi)^dat$n), 
                       dat$y*log(psi) + (dat$n - dat$y)*log(1-psi) + log(theta)))
  loglik_keep <- numeric(n.reps)
  loglik_keep[1] <- loglik
  
  #jumping std devs for metropolis (include these as an option?)
  jump_beta <- rep(0.10, n.players)
  jump_delta <- rep(0.10, n.players)
  jump_alpha <- rep(0.10, n.players)
  jump_lambda <- rep(0.10, n.players)
  A <- 1.1 
  B <- 1.1^(-44/56)
  for (k in 2:n.reps) {
    if (k%%100==0) cat(k, "\n")
    
    #metropolis step for \beta
    stuff <- metropolis(beta=beta, alpha=alpha, delta=delta, lambda=lambda,
                        jump_sd=jump_beta, param="beta")
    beta_keep[k,] <- beta <- stuff$beta
    if (tune & k < adapt) jump_beta <- ifelse(stuff$ratio > 0.44, jump_beta*A, jump_beta*B)
    
    #metropolis step for \delta
    stuff <- metropolis(beta=beta, alpha=alpha, delta=delta, lambda=lambda,
                        jump_sd=jump_delta, param="delta")
    delta_keep[k,] <- delta <- stuff$delta
    if (tune & k < adapt) jump_delta <- ifelse(stuff$ratio > 0.44, jump_delta*A, jump_delta*B)
    
    #metropolis step for \alpha
    stuff <- metropolis(beta=beta, alpha=alpha, delta=delta, lambda=lambda,
                        jump_sd=jump_alpha, param="alpha")
    alpha_keep[k,] <- alpha <- stuff$alpha
    if (tune & k < adapt) jump_alpha <- ifelse(stuff$ratio > 0.44, jump_alpha*A, jump_alpha*B)
    
    #metropolis step for \lambda
    stuff <- metropolis(beta=beta, alpha=alpha, delta=delta, lambda=lambda,
                        jump_sd=jump_lambda, param="lambda")
    lambda_keep[k,] <- lambda <- stuff$lambda
    if (tune & k < adapt) jump_lambda <- ifelse(stuff$ratio > 0.44, jump_lambda*A, jump_lambda*B)
    #calculate log-likelihood
    t <- exp(rep(beta, ni) + rep(alpha, ni)*dat$home)  #save repeated calc for theta
    theta <- t/(1 + t)
    p <- exp(rep(delta, ni) + rep(lambda, ni)*dat$home)  #save repeated calc for lambda
    psi <- p/(1 + p)
    loglik <- sum(ifelse(zero, log((1-theta) + theta*(1-psi)^dat$n), 
                         dat$y*log(psi) + (dat$n - dat$y)*log(1-psi) + log(theta)))
    loglik_keep[k] <- loglik
  }
  return(list(betas=beta_keep,
              deltas=delta_keep,
              alphas=alpha_keep,
              lambdas=lambda_keep, 
              logliks=loglik_keep))
}