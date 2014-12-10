# This script is a simple proof that the MCMC works
# That is, we can revover the true parameters from simulated data
source("code/01-dat.R")

#"known" hyperparameters
sigma_alpha.true <- 1  #eventually we'll learn this guy
sigma_lambda.true <- 1 #eventually we'll learn this guy
sigma_beta.true <- 1
sigma_delta.true <- 1
alpha0 <- 2
beta0 <- -2
delta0 <- -3
lambda0 <- 1

alphas.true <- rnorm(n.players, alpha0, 0.5)
betas.true <- rnorm(n.players, beta0, 0.5)
t <- exp(rep(betas.true, ni) + rep(alphas.true, ni)*dat$home)
theta.true <- t/(1+t)
x.sim <- rbinom(N, size=1, prob=theta.true)

lambdas.true <- rnorm(n.players, lambda0, 0.5)
deltas.true <- rnorm(n.players, delta0, 0.5)
p <- exp(rep(deltas.true, ni) + rep(lambdas.true, ni)*dat$home)
psi.true <- p/(1+p)
y.sim <- rbinom(N, size=dat$n, prob=psi.true)
y.sim[x.sim == 0] <- 0

loglike.true <- sum(ifelse(y.sim == 0, log((1-theta.true) + theta.true*(1-psi.true)^dat$n), 
                                       dat$y*log(psi.true) + (dat$n - dat$y)*log(1-psi.true) + log(theta.true)))

dat.sim <- dat
dat.sim$y <- y.sim

# Run the MCMC
source("code/02-mcmc.R")
n.iter <- 1000
sim.res <- mcmc_chain(dat.sim, reps=n.iter)







#recover betas
par(mfrow=c(2,2))
n.plots <- sample(1:n.players, 4)
result <- sim.res[[1]]
for (i in n.plots) {
  plot(1:n.iter, result$betas[,i], col=1, type="l")
  abline(h=betas.true[i], col=2, lty=2)
  if (length(sim.res) > 1) {
    for (j in 2:length(sim.res)) {
      result <- sim.res[[j]]
      lines(1:n.iter, result$betas[, i], col=j, type="l")
    }
  }
}

#recover alphas
par(mfrow=c(2,2))
n.plots <- sample(1:n.players, 4)
result <- sim.res[[1]]
for (i in n.plots) {
  plot(1:n.iter, result$alphas[,i], col=1, type="l")
  abline(h=alphas.true[i], col=2, lty=2)
  if (length(sim.res) > 1) {
    for (j in 2:length(sim.res)) {
      result <- sim.res[[j]]
      lines(1:n.iter, result$alphas[, i], col=j, type="l")
    }
  }
}

#recover deltas
par(mfrow=c(2,2))
n.plots <- sample(1:n.players, 4)
result <- sim.res[[1]]
for (i in n.plots) {
  plot(1:n.iter, result$deltas[,i], col=1, type="l")
  abline(h=deltas.true[i], col=2, lty=2)
  if (length(sim.res) > 1) {
    for (j in 2:length(sim.res)) {
      result <- sim.res[[j]]
      lines(1:n.iter, result$deltas[, i], col=j, type="l")
    }
  }
}

#recover lambdas
par(mfrow=c(2,2))
n.plots <- sample(1:n.players, 4)
result <- sim.res[[1]]
for (i in n.plots) {
  plot(1:n.iter, result$lambdas[,i], col=1, type="l")
  abline(h=lambdas.true[i], col=2, lty=2)
  if (length(res) > 1) {
    for (j in 2:length(res)) {
      result <- sim.res[[j]]
      lines(1:n.iter, result$lambdas[, i], col=j, type="l")
    }
  }
}
