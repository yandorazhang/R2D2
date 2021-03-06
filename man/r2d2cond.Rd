% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2d2cond.R
\name{r2d2cond}
\alias{r2d2cond}
\title{r2d2cond}
\usage{
r2d2cond(
  x,
  y,
  a,
  b,
  a1 = 0.01,
  b1 = 0.01,
  nu = NULL,
  mu = NULL,
  q = NULL,
  mcmc.n = 10000,
  gibbs = NULL
)
}
\arguments{
\item{x}{Scaled design matrix with no intercept (such that diag(x'x)=1).}

\item{y}{Centered response variable (mean 0).}

\item{a, b}{Shape parameters of Beta distribution for R-Squared.}

\item{a1, b1}{Shape and scale parameters of Inverse Gamma distribution on
sigma^2.}

\item{nu}{Shape parameter of Gamma distribution for local variance
components.}

\item{mu}{Rate parameter of Gamma distribution for local variance components
(optional).}

\item{q}{Number of expected non-zero coefficients (if nu was not specified
directly).}

\item{mcmc.n}{Number of MCMC iterations.}

\item{gibbs}{Boolean, whether to use Gibbs sampling. Applicable only for
uniform-on-ellipsoid prior (nu=mu=q=NULL), and a <= p/2. Will be used by
default if possible.}
}
\value{
A list with the following components:
\itemize{
\item{beta: Matrix (mcmc.n * p) of posterior samples.}
\item{sigma2: Vector (mcmc.n) of posterior samples.}
\item{theta: Vector (mcmc.n) of posterior samples, if applicable.}
\item{Lambda: Matrix (mcmc.n * p) of posterior samples, if applicable.}
\item{w.samples: Vector (mcmc.n) of posterior samples, if Gibbs sampler
was used.}
\item{z.samples: Vector (mcmc.n) of posterior samples, if Gibbs sampler
was used.}
\item{acc.rates: List of acceptance rates for sigma2, theta and beta.
}
}
}
\description{
\code{r2d2cond} adopts MCMC sampling algorithms, aiming to obtain a sequence
of posterior samples based on the conditional distribution on R-squared.
}
\section{Notes}{

\itemize{
\item{The uniform-on-ellipsoid prior will be used if nu, mu, and q are all
NULL.}
\item{To fit the local shrinkage model, nu or q must be specified, and mu
is optional (a suitable value will be chosen based on nu, p, and n).}
}
}

\examples{
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
}
