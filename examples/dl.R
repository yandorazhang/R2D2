rho <- 0.5
# Number of predictors
p <- 25
# Number of observations
n <- 60
# Construct beta
n_nonzero <- 5
beta <- rep(0, p)
beta[11:(10 + n_nonzero)] <- stats::rt(n_nonzero, df = 3) * sqrt(0.5/(3 * n_nonzero/2))
# Construct x
sigma <- 1
times <- 1:p
H <- abs(outer(times, times, "-"))
V <- sigma * rho^H
x <- mvtnorm::rmvnorm(n, rep(0, p), V)
x <- scale(x, center = TRUE, scale = FALSE)
# Construct y
y <- x %*% beta + stats::rnorm(n)

# Gibbs sampling
a.prior0 <- p/n
c.prior0 <- a.prior0/p
a.dl.prior0 <- 2 * c.prior0
mcmc.n <- 10000
fit.dl <- dl(x = x, y = y, hyper = list(a1 = 0.001, b1 = 0.001),
                  a.prior = a.dl.prior0, mcmc.n = mcmc.n, print = FALSE)

# Discard the early samples and get statistics of beta
burnIn <- 5000
beta.dl = fit.dl$beta[burnIn:mcmc.n, ]
mean.beta = apply(beta.dl, 2, mean)
