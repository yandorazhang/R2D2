#' Gibbs Sampler for Marginal R2-D2
#'
#' \code{r2d2marg} adopts the Gibbs sampling algorithm, aiming to obtain a
#' sequence of posterior samples, which are approximately from the Marginal
#' R2-D2 prior.
#'
#' @param x An n-by-p input matrix, each row of which is an observation vector.
#' @param y An n-by-one vector, representing the response variable.
#' @param hyper The values of the hyperparameters in the prior, i.e., the values
#'   of \code{(a_pi, b, a1, b1)}.
#'   \itemize{
#'   \item{\code{a_pi} is the concentration parameter of the Dirichlet
#'   distribution, which controls the local shrinkage of phi. Smaller values of
#'   a_pi lead to most phi close to zero; while larger values of a_pi lead to a
#'   more uniform phi.}
#'   \item{\code{b} is the second shape parameter of beta prior for the
#'   R-squared.}
#'   \item{\code{a1} and \code{b1} are shape and scale parameters of Inverse
#'   Gamma distribution on sigma^2.} }
#'   \code{hyper} is set to be (1/(p^(b/2)n^(b/2)logn), 0.5, 0.001, 0.001) by
#'   default.
#' @param mcmc.n Number of MCMC samples, by default, 10000.
#' @param eps Tolerance of convergence, by default, 1e-7.
#' @param thin Thinning parameter of the chain. The default is 1, representing
#'   no thinning.
#' @param print Boolean variable determining whether to print the progress of
#'   the MCMC iterations, i.e., which iteration the function is currently on.
#'   The default is TRUE.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item{beta: Matrix (mcmc.n * p) of posterior samples for beta.}
#'   \item{sigma2: Vector (mcmc.n) of posterior samples for sigma^2.}
#'   \item{psi: Matrix (mcmc.n * p) of posterior samples for psi.}
#'   \item{w: Vector (mcmc.n) of posterior samples for the total prior
#'   probability w.}
#'   \item{xi: Vector (mcmc.n) of posterior samples for xi.}
#' }
#'
#' @example inst/examples/r2d2marg.R
#'
#' @export
r2d2marg <- function(x, y, hyper, mcmc.n = 10000, eps = 1e-07, thin = 1, print = TRUE) {

    #------------------------------------------
    EPS = 1e+20
    # 1. n, p, max(n,p)
    p <- ncol(x)
    n <- nrow(x)
    max.np <- max(n, p)

    # 2. hyperparameter values (a_pi, b, a1, b1)
    if (missing(hyper)) {
        b <- 0.5
        a_pi <- 1/(p^(b/2) * n^(b/2) * log(n))
        a1 <- 0.001
        b1 <- 0.001
    } else {
        a_pi <- hyper$a_pi
        b <- hyper$b
        a1 <- hyper$a1
        b1 <- hyper$b1
    }

    a <- a_pi * p

    # 3. Define variables to store posterior samples
    keep.beta = keep.phi = keep.psi = matrix(0, nrow = mcmc.n, ncol = p)
    keep.w = keep.xi = keep.sigma2 = rep(0, mcmc.n)

    # 4. Initial values.
    ## beta
    tem <- stats::coef(stats::lm(y ~ x - 1))
    tem[which(is.na(tem))] <- eps
    beta <- tem
    ## sigma
    sigma2 <- MCMCpack::rinvgamma(1, shape = a1 + n/2, scale = b1 + 0.5 * crossprod(y - x %*% beta))
    sigma <- sqrt(sigma2)
    ## phi
    phi <- rep(a_pi, p)
    ## xi
    xi <- a/(2 * b)
    ## w
    w <- b
    ## psi
    psi <- stats::rexp(p, rate = 0.5)
    #------------------------------------------

    XTX <- crossprod(x)  #t(X)%*%X
    XTY <- crossprod(x, y)  #t(X)%*%y

    keep.beta[1, ] <- beta
    keep.sigma2[1] <- sigma2
    keep.phi[1, ] <- phi
    keep.w[1] <- w
    keep.xi[1] <- xi
    keep.psi[1, ] <- psi

    #------------------------------------------
    #------------------------------------------
    # Now MCMC runnings!
    for (i in 2:mcmc.n) {
        for (j in 1:thin) {
            a <- a_pi * p

            # draw a number from N(0,1)
            Z <- stats::rnorm(p, mean = 0, sd = 1)

            #------------------------------------------
            # (i) Sample beta| phi, w, sigma2, y
            d.inv <- 1/(psi * phi * w)
            ad.inv <- abs(d.inv)
            sd.inv <- sign(d.inv)

            inx.e <- which(ad.inv < eps)
            inx.E <- which(ad.inv > EPS)
            d.inv[inx.e] <- eps * sd.inv[inx.e]
            d.inv[which(is.infinite((d.inv)))] <- EPS
            d.inv[inx.E] <- EPS * sd.inv[inx.E]

            if (length(d.inv) == 1) {
                Dinv <- d.inv
            } else {
                Dinv <- diag(d.inv)
            }

            Vinv <- XTX + Dinv

            # ********************************#********************************
            temQ <- chol(Vinv, pivot = T, tol = 1e-100000)
            pivot <- attr(temQ, "pivot")
            temc <- temQ[, order(pivot)]
            B <- XTY + t(temc) %*% Z * sigma


            # s = svd(Vinv) tes = s$u %*% diag(sqrt(s$d)) %*% t(s$v) b = XTY + t(tes)%*%Z*sqrt(sigma2[k-1])

            # b = XTY + t(chol(Vinv))%*%Z*sqrt(sigma2[k-1])
            # ********************************#********************************

            beta <- solve(Vinv, B, tol = 1e-10000)
            beta[which(is.na(beta))] <- 0

            #------------------------------------------
            # (ii) Sample sigma2| beta, phi, w, y
            sigma2 <- MCMCpack::rinvgamma(1, shape = a1 + n/2 + p/2, scale = b1 + 0.5 * crossprod(y - x %*%
                beta) + sum(beta * d.inv * beta)/2)
            sigma <- sqrt(sigma2)

            #------------------------------------------
            # (iii) Sample psi| beta, phi, w, sigma2
            #------------------------------------------
            mu.te <- sigma * sqrt(phi * w)/abs(beta)
            amu.te <- abs(mu.te)
            smu.te <- sign(mu.te)
            inx.e <- which(amu.te < eps)
            inx.E <- which(amu.te > EPS)
            mu.te[inx.E] <- EPS * smu.te[inx.E]
            mu.te[inx.e] <- eps * smu.te[inx.e]
            ## rinvgauss_new: self-defined function
            psi <- 1/rinvgauss_new(p, mean = mu.te, shape = 1)


            #------------------------------------------
            # (iv) Sample w|beta, phi, xi, sigma2
            chi.te <- sum(beta^2/(psi * phi))/sigma2
            achi.te <- abs(chi.te)
            schi.te <- sign(chi.te)
            inx.e <- which(achi.te < eps)
            inx.E <- which(achi.te > EPS)
            chi.te[inx.e] <- eps * schi.te[inx.e]
            chi.te[inx.E] <- EPS * schi.te[inx.E]
            chi.te[which(chi.te == 0)] <- eps
            w <- GIGrvg::rgig(n = 1, lambda = a - p/2, chi = chi.te, psi = 4 * xi)



            #------------------------------------------
            # (v) Sample xi|w
            xi <- stats::rgamma(1, shape = a + b, rate = 1 + 2 * w)



            #------------------------------------------
            # (vi) Sample phi|beta, xi, sigma2, y
            TT <- rep(0, p)
            tem <- (beta^2/sigma2)/psi
            atem <- abs(tem)
            stem <- sign(tem)
            inx.e <- which(atem < eps)
            inx.E <- which(atem > EPS)
            tem[inx.e] <- eps * stem[inx.e]
            tem[inx.E] <- EPS * stem[inx.E]
            tem[which(tem == 0)] <- eps
            TT <- apply(matrix(tem, ncol = 1), 1, function(xx) return(GIGrvg::rgig(n = 1, lambda = a_pi -
                0.5, chi = xx, psi = 4 * xi)))

            Ts <- sum(TT)
            phi <- TT/Ts

        }

        if (print) {
            if (i%%500 == 0) {
                print(paste(c("The ", i, "th sample."), collapse = ""))
            }
        }

        keep.beta[i, ] <- beta
        keep.sigma2[i] <- sigma2
        keep.phi[i, ] <- phi
        keep.w[i] <- w
        keep.xi[i] <- xi
        keep.psi[i, ] <- psi

        if (sqrt(sum((keep.beta[i, ] - keep.beta[i - 1, ])^2)) < eps & abs(keep.xi[i] - keep.xi[i - 1]) <
            eps & sqrt(sum((keep.phi[i, ] - keep.phi[i - 1, ])^2)) < eps & abs(keep.w[i] - keep.w[i -
            1]) < eps & abs(keep.sigma2[i] - keep.sigma2[i - 1]) < eps & sqrt(sum((keep.psi[i, ] - keep.psi[i -
            1, ])^2)) < eps)
            break

    }

    return(list(beta = keep.beta, sigma2 = keep.sigma2, psi = keep.psi, w = keep.w, xi = keep.xi))
}

