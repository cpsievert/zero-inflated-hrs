library("coda")
library("xtable")

# Make sure we have the data and the run_mcmc function
source("code/01-dat.R")
source("code/02-mcmc.R")

# Set some MCMC options
n.iter <- 10000
adapt <- 500
burnin <- 5000

res <- mcmc_chain(dat, reps=n.iter)
# This results object is used to create figures in the paper
#save(res, file="results.rda")

alphas <- mcmc.list(lapply(res, function(x) mcmc(x$alphas[burnin:n.iter,])))
betas <- mcmc.list(lapply(res, function(x) mcmc(x$betas[burnin:n.iter,])))
deltas <- mcmc.list(lapply(res, function(x) mcmc(x$deltas[burnin:n.iter,])))
lambdas <- mcmc.list(lapply(res, function(x) mcmc(x$lambdas[burnin:n.iter,])))
# convergence diagnostics -- should be close to 1
gelman.diag(alphas)
gelman.diag(betas)
gelman.diag(deltas)
gelman.diag(lambdas)

acfplot(alphas)
acfplot(betas)
acfplot(deltas) #yikes
acfplot(lambdas) #yikes

t <- exp(as.matrix(betas))
theta <- t/(1+t)
hist(theta[,1])  #wtf

p <- exp(as.matrix(lambdas))
psi <- p/(1+p)
hist(psi[,1])   #wtf

#95% credible intervals
summary(alphas)[[2]][,c(1,5)]
summary(lambdas)[[2]][,c(1,5)]

betaz <- apply(as.matrix(betas), 1, function(x) rep(x, ni))
alphaz <- apply(as.matrix(alphas), 1, function(x) rep(x, ni))
t <- exp(betaz + alphaz*dat$home)
thetaz <- t/(1+t)

lambdaz <- apply(as.matrix(lambdas), 1, function(x) rep(x, ni))
deltaz <- apply(as.matrix(deltas), 1, function(x) rep(x, ni))
p <- exp(deltaz + lambdaz*dat$home)
psiz <- p/(1+p)

# we *really* want the posterior distribution of *home runs* given potential
# then, ecdf tells us the probability of hitting x or less home runs

#cumulative dist (away)
par(mfrow=c(1,1))
edwina <- ecdf(theta[,4])
plot(edwina, main="", ylab="Empirical Distribution Function", 
     xlab="Home run potential")
edwinh <- ecdf(thetaz[ni[4]+1,])
lines(edwinh, col=2)
legend(0.8, 0.4, c("Home", "Away"), col=2:1, lty=1)

#cumulative dist (away)
par(mfrow=c(1,1))
josha <- ecdf(4 *psi[,2])
plot(josha, main="", ylab="Empirical Distribution Function", 
     xlab="Home runs given potential in a game with 4 atbats")
joshh <- ecdf(4 * psiz[ni[2]+5,])
lines(joshh, col=2)
legend(3, 0.4, c("Home", "Away"), col=2:1, lty=1)


#posterior p-values
#pg 312 of the notes suggests checking the "variance or perhaps the range"
n.zeros <- NULL
ranges <- NULL
for (m in 1:1000) {  
  sim.x <- rbinom(N, size=1, prob=theta[m,])
  sim.y <- rbinom(N, size=dat$n, prob=psi[m,])
  sim.y[sim.x == 0] <- 0
  n.zeros <- c(n.zeros, sum(sim.y == 0))
  ranges <- c(ranges, max(sim.y) - min(sim.y))
}

real.n.zero <- sum(dat$y == 0)
p1 <- sum(n.zeros > real.n.zero)/length(n.zeros)
real.range <- max(dat$y)-min(dat$y)
p2 <- sum(ranges > real.range)/length(ranges)


hist(n.zeros)
abline(v=real.n.zero, col=2)
hist(ranges)

abline(v=real.range, col=2)


#acceptance probs
lu <- function(x) length(unique(x))
ac.alpha <- apply(res[[1]]$alphas, 2, lu)/n.iter
ac.beta <- apply(res[[1]]$betas, 2, lu)/n.iter
ac.delta <- apply(res[[1]]$deltas, 2, lu)/n.iter
ac.lambda <- apply(res[[1]]$lambda, 2, lu)/n.iter
#hist(ac.alpha)
#hist(ac.beta)
#hist(ac.delta)
#hist(ac.lambda)

df <- data.frame(alphas=ac.alpha, betas=ac.beta, deltas=ac.delta, lambdas=ac.lambda)
xtable(df, caption="Acceptance Rates for model parameters sampled via Metropolis-Hastings")
