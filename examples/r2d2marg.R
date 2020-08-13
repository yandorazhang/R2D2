rho <- 0.5
# Number of predictors
p <- 25
# Number of observations
n <- 60
# Construct beta
n_nonzero <- 5
beta <- rep(0, p)
set.seed(1)
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
mcmc.n <- 10000
fit.new <- r2d2marg(x = x, y = y, mcmc.n = mcmc.n, print = FALSE)

# Discard the early samples
burnIn <- 5000
beta.new <- fit.new$beta[burnIn:10000, ]
