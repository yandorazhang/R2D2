rMvNormLeft = function(mu, SigmaHalf) {
  # Sample from a MVN(mu, SigmaHalf %*% t(SigmaHalf)) distribution. Faster than
  # MASS::mvrnorm, if the sqrt matrix is already available.
  # Args:
  #   mu: mean vector.
  #   SigmaHalf: left square root matrix of Sigma
  return(as.vector(mu + SigmaHalf %*% stats::rnorm(dim(SigmaHalf)[2])))
}

rMvNormChol = function(mu, Sigma) {
  # Sample from a MVN(mu, Sigma) distribution using cholesky decomposition
  # with pivoting. Works when MASS::mvrnorm doesn't.
  # Args:
  #   mu: mean vector
  #   Sigma: covariance matrix
  R.pivoted <- suppressWarnings(chol(Sigma, pivot = T)) # warning if not PD
  oo <- order(attr(R.pivoted, "pivot"))
  SigmaHalf <- t(R.pivoted[,oo])
  return(rMvNormLeft(mu, SigmaHalf))
}

rBingham = function(n, A) {
  # Borrowed heavily from Directional::rbingham(), but uses optimal bound M(b)
  # from Kent, Ganeiber and Mardia (2013), and some speed ups to sampling if A
  # is a vector (diag(A)) then skip eigen, but assume it is sorted (decr)
  if (is.vector(A)) {
    p <- length(A)
    lam <- A
    V <- diag(p)
  } else {
    p <- ncol(A)
    eig <- eigen(A, symmetric = T)
    lam <- eig$values
    V <- eig$vectors
  }
  lam <- lam - lam[p]
  lam <- lam[-p]
  # lam <- sort(lam, decreasing = TRUE) # already sorted
  nsamp <- 0
  X <- NULL
  lam.full <- c(lam, 0)
  qa <- length(lam.full)
  mu <- numeric(qa)

  FunB = function(b,lam) {
    # funtion of b that is solution to minimizing M(b)
    sum(1/(b+2*lam))
  }
  FindB = function(lam.full, tol=1e-4) {
    # Uses bisection search to find the solution to FunB(b) = 1
    max.b = 1
    while (FunB(max.b,lam.full)>1) {
      max.b = max.b*2
    }
    min.b = 0
    b = (min.b + max.b)/2
    target = FunB(b,lam.full) - 1
    while(abs(target) > tol) {
      if (target > 0) {
        min.b = b
      } else {
        max.b = b
      }
      b = (min.b + max.b)/2
      target = FunB(b,lam.full) - 1
    }
    return(b)
  }
  b = FindB(lam.full)
  # note: b=1 in Directional::rbingham()

  sigacginv <- 1 + 2 * lam.full / b
  Ntry <- 0
  while (nsamp < n) {
    x.samp <- FALSE
    while (x.samp == FALSE) {
      yp <- stats::rnorm(n = qa, mean=0, sd = 1/sqrt(sigacginv))
      y <- yp/sqrt(sum(yp^2))
      lratio <- -sum(y^2 * lam.full) + qa/2 * log(b/qa) +
        0.5 * (qa - b) + qa/2 * log(1+2/b*sum(y^2 * lam.full))
      if (log(stats::runif(1)) < lratio) {
        X <- c(X, y)
        x.samp <- TRUE
        nsamp <- nsamp + 1
      }
      Ntry <- Ntry + 1
    }
  }
  x <- matrix(X, byrow = TRUE, ncol = qa)
  tcrossprod(x, V)
}

SampleFB = function(kappa.mu, A, max.iters=5000) {
  # Rejection sampler to sample from Fisher-Bingham(kappa*mu, A) based on
  # Kent, et al (2013) "A new method to simulate the Bingham..."
  # f(x) = C * exp(kappa.mu'x - x'Ax)
  # NOTE: Directional::rfb is only for p=3 and is wrong
  # Args:
  #   kappa.mu: vector
  #   A: matrix
  #   max.iters: maximum number of iterations for rejection sampler.
  p <- length(kappa.mu)
  k <- sqrt(sum(kappa.mu^2)) # kappa
  m <- kappa.mu / k          # mu
  A.1 <- k/2*(diag(p) - tcrossprod(m)) + A
  i <- 1
  # Take eigendecomp once outside of loop:
  e.A1 <- eigen(A.1, symmetric = T)
  A.1V <- e.A1$vectors
  A.1d <- e.A1$values
  while (i <= max.iters) {
    x <- as.vector(tcrossprod(rBingham(1, A.1d), A.1V))
    a <- exp(-(1-sum(m*x))^2 * k/2)
    if (stats::runif(1) < a) {
      return (list(x=x, i=i, finish=T))
    }
    i <- i+1
  }
  return (list(x=x, i=i, finish=F))
}

CalcLogPost <- function(X, y, Beta, sigma2, n, p, a, b, alpha.0, beta.0) {
  # Calculates log-posterior for Metropolis-Hastings updating for (non-local
  # shrinkage model).
  XBeta = X %*% Beta
  BXXBn = sum(XBeta^2) / n
  ymXBeta = y - XBeta
  sse = sum(ymXBeta^2)
  log.post = (a-p/2) * log(BXXBn) - (a+b) * log(1+BXXBn/sigma2)
  log.post = log.post - (n/2+a+alpha.0+1) * log(sigma2)
  log.post = log.post - (sse/2 + beta.0) / sigma2
  return (log.post)
}

CalcLogPostS2Th <- function(ypy, ypX, VDhalfinv, Gam, sigma2, theta, a, b,
                            alpha.0, beta.0, n, p) {
  # Calculates log-posterior for sigma2 and theta for M-H updating.
  log.post <- (a-1) * log(theta) - (a+b) * log(1+theta)
  log.post <- log.post - (alpha.0+1+n/2) * log(sigma2)
  log.post <- log.post - (ypy/2 + beta.0) / sigma2 - n*theta/2
  log.post <- log.post + as.vector(ypX %*%(VDhalfinv%*%Gam)*sqrt(theta/sigma2))
  return (log.post)
}

TuneBeta <- function(Beta.var, tmpacc.beta, tmpatt) {
  if (tmpacc.beta / tmpatt < .1)  Beta.var <- .8 * Beta.var
  if (tmpacc.beta / tmpatt > .45) Beta.var <- 1.2 * Beta.var
  return (Beta.var)
}

TuneSigma <- function(sigma2.var, tmpacc.s2, tmpatt) {
  if (tmpacc.s2 / tmpatt < .35) sigma2.var <- .8 * sigma2.var
  if (tmpacc.s2 / tmpatt > .55) sigma2.var <- 1.2 * sigma2.var
  return (sigma2.var)
}

GetMu <- function(nu, n, p, alpha = .01, init.mu = 1, max.iters = 1000,
                  n.samples = 1e6) {
  # Bisection search to find mu given nu, n and p.
  tol <- min(.001, alpha / 10)
  Y <- stats::rgamma(n.samples, nu, 1)
  vals <- 1 / (2 * Y)
  mu <- init.mu
  low <- 0
  for (i in 1:max.iters) {
    obj <- mean(stats::pgamma(vals, 0.5, mu)) - (1 - alpha)
    if (abs(obj) < tol) {
      break
    } else {
      if (obj > 0) {
        # mu too big
        mu <- (low + mu) / 2
      } else {
        # mu too small
        low <- mu
        mu <- 2 * mu
      }
    }
  }
  return (mu)
}

GetNu <- function(p, q, C = 0.95, v=1, nu.init = 1, N = 1e4, eps = 0.001) {
  # Finds nu, such that the median sum of the largest q order statistics of
  #   a sample from a draw from independent gamma(nu, v) r.v.'s is
  #   approximately equal to C.
  #   Note: v doesn't matter
  if (q/p > C) stop(paste(round(q/p,2),'must be less than',C))
  nu = nu.init
  iters <- 1
  ub <- NA
  lb <- 0
  repeat {
    X <- matrix(stats::rgamma(N*p, nu, v), nrow=N, ncol=p)
    totals <- rowSums(X)
    if (q == 1) {
      partial.sums <- apply(X, 1, max)
    } else {
      partial.sums = colSums(apply(X, 1, sort, decreasing = T)[1:q,])
    }
    C.hat <- stats::median(partial.sums / totals)
    if (abs(C.hat - C) < eps) break
    if (C.hat < C) {
      ub <- nu
      nu <- mean(c(lb, ub))
    } else {
      lb <- nu
      if (is.na(ub)) {
        nu <- nu * 2
      } else {
        nu <- mean(c(lb, ub))
      }
    }
    iters <- iters + 1
  }
  return (nu)
}
