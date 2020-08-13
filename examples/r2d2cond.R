# Load cereal data and standardize n = 15, p =145
data("cereal", package = "chemometrics")
X <- scale(cereal$X)/sqrt(nrow(cereal$X) - 1)
y <- scale(cereal$Y[, "Starch"], scale = FALSE)
n <- nrow(X)
p <- ncol(X)

# Number of iterations and burn-in period
burn.in <- 500
mcmc.n <- 10000

# Fit various priors and compute posterior means of Beta. Note that since p > n
# for this data, the Uniform-on-ellipsoid is not valid. We demonstrate the
# approach using only the first 5 predictors here.

# Uniform-on-ellipsoid prior with a = b = 1.
fit.uoe <- r2d2cond(x = X[, 1:5], y = y, a = 1, b = 1, mcmc.n = mcmc.n)
beta.uoe <- colMeans(fit.uoe$beta[burn.in:mcmc.n, ])

# Uniform-on-ellipsoid prior with a = p/n and b = 0.5.
fit.uoe1 <- r2d2cond(x = X[, 1:5], y = y, a = 5/n, b = 0.5, mcmc.n = mcmc.n)
beta.uoe1 <- colMeans(fit.uoe1$beta[burn.in:mcmc.n, ])

# Compare with OLS
beta.ols <- coef(lm(y ~ X[, 1:5] - 1))

# print out coefficient estimates from 3 methods.
beta.uoe
beta.uoe1
beta.ols

# Now use full p = 145 with Local Shrinkage Model (as done in real data example
# with cereal data)

# Local shrinkage prior with a = b = 1.
fit.shrink <- r2d2cond(x = X, y = y, a = 1, b = 1, q = nrow(X), mcmc.n = mcmc.n)
beta.shrink <- colMeans(fit.shrink$beta[burn.in:mcmc.n, ])

# Local shrinkage prior with a = p/n, b = 0.5.
fit.shrink1 <- r2d2cond(x = X, y = y, a = p/n, b = 0.5, q = nrow(X), mcmc.n = mcmc.n)
beta.shrink1 <- colMeans(fit.shrink1$beta[burn.in:mcmc.n, ])
