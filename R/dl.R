#' Gibbs Sampler for Dirichlet-Laplace Prior
#'
#' \code{dl} aims to generate the posterior samples for Dirichlet-Laplace
#' prior.
#'
#' @param x An n-by-p input matrix, each row of which is an observation vector.
#' @param y An n-by-one vector, representing the response variable.
#' @param hyper The values of hyperparameters \code{(a1,b1)} in the prior of
#'   sigma^2.
#' @param a.prior The value of hyperparameter in the prior, which can be a value
#'   between (1/max(n,p), 1/2). By default, it is set at 1/max(n,p). 
#' @param mcmc.n Number of MCMC samples, by default, 10000.
#' @param thin Thinning parameter of the chain. The default is 1, representing
#'   no thinning.
#' @param eps Tolerance of convergence, by default, 1e-7.
#' @param print Boolean variable determining whether to print the progress of
#'   the MCMC iterations, i.e., which iteration the function is currently on.
#'   The default is TRUE.
#'
#' @return A list containing the following components:
#'   \itemize{
#'   \item{beta: Matrix (mcmc.n * p) of posterior samples for beta.}
#'   \item{psi: Matrix (mcmc.n * p) of posterior samples for psi.}
#'   \item{phi: Matrix (mcmc.n * p) of posterior samples for phi.}
#'   \item{tau: Vector (mcmc.n) of posterior samples for tau.}
#'   \item{sigma2: Vector (mcmc.n) of posterior samples for sigma^2.}
#'   }
#'
#' @example inst/examples/dl.R
#'
#' @export
dl <- function(x, y, hyper, a.prior, mcmc.n = 10000, thin = 1, eps = 1e-07, print = TRUE) {

    EPS = 1e+20
    # n, p, max(n,p)
    p <- ncol(x)
    n <- nrow(x)
    max.np <- max(n, p)

    # priors IG(a1,b1) for sigma2.
    if (missing(hyper)) {
        a1 <- 0.001
        b1 <- 0.001
    } else {
        a1 <- hyper$a1
        b1 <- hyper$b1
    }

    # the discrete uniform value for a.
    if (missing(a.prior)) {
        a.prior <-  1/max.np 
    }


    # Define variables to store posterior samples
    beta = phi = psi = matrix(0, nrow = mcmc.n, ncol = p)
    tau = sigma2 = rep(0, mcmc.n)
    #------------------------------------------
    #------------------------------------------
    # Initial values.
    tem <- stats::coef(stats::lm(y ~ x - 1))
    tem[which(is.na(tem))] <- eps
    beta[1, ] <- tem
    sigma2[1] <- MCMCpack::rinvgamma(1, shape = a1 + n/2, scale = b1 + 0.5 * crossprod(y - x %*% beta[1, ]))

    phi[1, ] <- rep(1/p, p)
    psi[1, ] <- stats::rexp(p, 0.5)
    tau[1] <- stats::rgamma(1, shape = p * a.prior, rate = 0.5)
    

    #------------------------------------------
    #------------------------------------------

    XTX <- crossprod(x)  #t(X)%*%X
    XTY <- crossprod(x, y)  #t(X)%*%y


    # Z <- rmvnorm(mcmc.n,mean=rep(0,p),sigma=diag(p))

    # Now MCMC runnings!
    for (k in 2:mcmc.n) {
        for (j in 1:thin) {

            Z <- stats::rnorm(p, mean = 0, sd = 1)

            #------------------------------------------
            # (ii) Sample beta| psi, phi, tau, y, sigma2
            d.te <- 1/(psi[k - 1, ] * (phi[k - 1, ]^2) * (tau[k - 1]^2))
            ad.te <- abs(d.te)
            sd.te <- sign(d.te)
            inx.e <- which(ad.te < eps)
            inx.E <- which(ad.te > EPS)
            d.te[inx.e] <- eps * sd.te[inx.e]
            d.te[which(is.infinite((d.te)))] <- EPS
            d.te[inx.E] <- EPS * sd.te[inx.E]


            if (length(d.te) == 1) {
                Dinv <- d.te
            } else {
                Dinv <- diag(d.te)
            }

            Vinv <- XTX + Dinv
            A <- Vinv

            # ********************************#********************************
            temQ <- chol(Vinv, pivot = T, tol = 0)
            pivot <- attr(temQ, "pivot")
            temc <- temQ[, order(pivot)]
            b <- XTY + t(temc) %*% Z * sqrt(sigma2[k - 1])  #*******


            # s = svd(Vinv) tes = s$u %*% diag(sqrt(s$d)) %*% t(s$v) b = XTY + t(tes)%*%Z*sqrt(sigma2[k-1]) #*******
            # ********************************#********************************


            # b = XTY + t(chol(Vinv))%*%Z*sqrt(sigma2[k-1]) #*******

            beta[k, ] <- solve(A, b, tol = 0)  #*******
            beta[k, which(is.na(beta[k, ]))] <- 0
            # print(system.time(chol(Vinv))) print(system.time(solve(A,b,tol=1e-1000)))


            #------------------------------------------
            # (i) Sample sigma2| beta, y
            sigma2[k] <- MCMCpack::rinvgamma(1, shape = a1 + n/2 + p/2, scale = b1 + 0.5 * crossprod(y - x %*% beta[k, ]) + sum(beta[k,
                ] * d.te * beta[k, ])/2)  #*******


            sigma.k <- sqrt(sigma2[k])

            #------------------------------------------
            # (iii) Sample psi|phi,tau,beta
            mu.te <- (phi[k - 1, ]/abs(beta[k, ]) * sigma.k) * tau[k - 1]  #*******
            amu.te <- abs(mu.te)
            smu.te <- sign(mu.te)
            inx.e <- which(amu.te < eps)
            inx.E <- which(amu.te > EPS)
            mu.te[inx.E] <- EPS * smu.te[inx.E]
            mu.te[inx.e] <- eps * smu.te[inx.e]

            psi[k, ] <- 1/rinvgauss_new(p, mean = mu.te, shape = 1)  #*******

            #------------------------------------------
            # (iv) Sample tau|phi,beta
            chi.te <- 2 * sum(abs(beta[k, ])/phi[k - 1, ])/sigma.k  #*******
            achi.te <- abs(chi.te)
            schi.te <- sign(chi.te)
            inx.e <- which(achi.te < eps)
            inx.E <- which(achi.te > EPS)
            chi.te[inx.e] <- eps * schi.te[inx.e]
            chi.te[inx.E] <- EPS * schi.te[inx.E]
            chi.te[which(chi.te == 0)] <- eps
            tau[k] <- rgig1(n = 1, lambda = p * a.prior - p, chi = chi.te, psi = 1)

            #------------------------------------------
            # (v) Sample phi|beta
            TT <- rep(0, p)
            tem <- 2 * abs(beta[k, ])/sigma.k  #*******
            atem <- abs(tem)
            stem <- sign(tem)
            inx.e <- which(atem < eps)
            inx.E <- which(atem > EPS)
            tem[inx.e] <- eps * stem[inx.e]
            tem[inx.E] <- EPS * stem[inx.E]
            tem[which(tem == 0)] <- eps
            TT <- apply(matrix(tem, ncol = 1), 1, function(xx) return(rgig1(n = 1, lambda = a.prior - 1, chi = xx, psi = 1)))

            # TT = apply(matrix(1:p,ncol=1),1,function(j) return(rgig(1,chi=2*tem[j], psi=1, lambda=a[k-1]-1)) ) te = matrix(2*tem,
            # ncol=1)

            # print(system.time(apply(matrix(1:p,ncol=1),1,function(j) return(rgig(1,chi=2*tem[j], psi=1, lambda=a[k-1]-1)) )) )
            # print(system.time(apply(te,1, function(xx) return(rgig(1, chi=xx, psi=1, lambda=a[k-1]-1)) )))

            Ts <- sum(TT)
            phi[k, ] <- TT/Ts

        }

        # #Plot the results thus far: if(k%%500==0){ plot(phi[1:k,1],type='l') }

        if (print) {
            if (k%%500 == 0) {
                print(paste(c("The ", k, "th sample."), collapse = ""))
            }
        }

        if (sqrt(sum((beta[k, ] - beta[k - 1, ])^2)) < eps & sqrt(sum((psi[k, ] - psi[k - 1, ])^2)) < eps & sqrt(sum((phi[k, ] -
            phi[k - 1, ])^2)) < eps & abs(tau[k] - tau[k - 1]) < eps & abs(sigma2[k] - sigma2[k - 1]) < eps)
            break
    }

    return(list(beta = beta[(1):mcmc.n, ], psi = psi[(1):mcmc.n, ], phi = phi[(1):mcmc.n, ], tau = tau[(1):mcmc.n], sigma2 = sigma2[(1):mcmc.n]))
}


rgig1 <- function(n = 1, lambda, chi, psi) {
    n = as.integer(n)
    lambda = as.double(lambda)
    chi = as.double(chi)
    psi = as.double(psi)
    GIGrvg::rgig(n, lambda, chi, psi)
}
