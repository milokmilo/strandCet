#' @title
#' Heligman-Pollard parameter coversion to age-specific probabilites of death.
#'
#' @description
#' Calculates the age-specific probabilities of death using the Heligman-Pollard model.
#'
#' @param theta A vector containing values for the 9 parameters of the adapted Heligman-Pollard model.
#' @param x A vector containing the ages at which to calculate the probabilities of death.
#'
#' @return
#' A vector of probabilities of death at ages defined by x.
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' @keywords Heligman-Pollard mortality bycatch
#'
#' @examples
#'
#' lifeN <- life.tab(cetaceans)
#'
#' modSi <- Si.mod(data = cetaceans, rm = 2,
#'                 par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))
#'
#' dataSi <- Si.pred(data = cetaceans, Sout = modSi, rm = 2)
#'
#' priors <- data.frame(priors.lo = c(0,0.5,0,0,0,0,6,0,1),
#'                      priors.hi = c(0.1,1,1,0.15,0.15,50,10,0.01,1.5))
#'
#' q0 <- HP.priors(pri.lo = priors$priors.lo,
#'                 pri.hi = priors$priors.hi,
#'                 theta.dim = 9)
#'
#' HP.mod(prior = q0, lifeTab = lifeN, rm = 2, K = 10, d = 10, B = 10, CI = 90)
#'
#' @export

mod <- function (theta, x) {
        A <- theta[1]
        B <- theta[2]
        C <- theta[3]
        D2 <- theta[4]
        D <- theta[5]
        E <- theta[6]
        F <- theta[7]
        G <- theta[8]
        H <- theta[9]
        f.x <- A^((x + B)^C) + (D2 + D * exp(-E * (log(x) - log(F))^2)) +
                (G * (H^x))/(1 + G * (H^x))
        f.x2 <- f.x
        f.x2[1e-06 > f.x] <- 1e-06
        f.x2[0.999999 < f.x] <- 0.999999
        return(f.x2)
}



#' @title
#' Binomial likelihood.
#'
#' @description
#' Calculates a log likelihood from a binomial distribution.
#'
#' @param x Same as x from \link{dbinom} - successes.
#' @param n Same as size from \link{dbinom} - trials.
#' @param p Same as prob from \link{dbinom} - observed probability.
#'
#' @references
#' R Development Core Team. (2009). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing Vienna, Austria.
#'
#' @seealso \link{dbinom}
#'
#' @keywords Heligman-Pollard binomial mortality bycatch
#'
#' @importFrom stats dbinom
#' @export

ll.binom <- function(x, n, p) {
        sum(dbinom(x = x, size = n, prob = p, log = TRUE))
}



#' @title
#' Prior likelihoods and weights.
#'
#' @description
#' Calculates the log-likelihood and importance weight for each set (i.e. each row) of Heligman-Pollard parameters in the prior.
#'
#' @param prior A ((theta.dim * 1000) * theta.dim) matrix containing the prior distribution.
#' @param nrisk A vector containing the number of persons at risk in each age group.
#' @param ndeath A vector containing the number of deaths in each age group.
#' @param theta.dim Number of columns of the prior matrix.
#' @param age A vector containing the ages at which each age interval begins.
#'
#' @return
#' \item{wts.0}{A vector containing an importance weight for each set of parameters from the prior.}
#' \item{log.like.0}{A vector containing a log likelihood for each set of parameters from the prior.}
#'
#' @note
#' Used in the \link{loop.optim} function
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' Poole, D and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Gen- eralized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{mod} \link{ll.binom} \link{loop.optim} \link{HP.mod}
#'
#' @keywords Heligman-Pollard priors mortality bycatch
#' @export

prior.likewts <- function (prior, nrisk, ndeath, age,
                           theta.dim = 9) {
        lx <- nrisk
        dx <- ndeath
        mp.mle <- function(theta, x.fit = age) {
                p.hat <- mod(theta = theta, x = x.fit)
                ll <- ll.binom(x = ndeath, n = nrisk, p = p.hat)
                return(ll)
        }
        B0 <- theta.dim * 1000
        log.like.0 <- rep(NA, B0)
        for (i in 1:B0) {
                log.like.0[i] <- mp.mle(prior[i, ])
        }
        a0 <- -max(log.like.0, na.rm = TRUE)
        like.0 <- exp(log.like.0 + a0 + 700)
        wts.0 <- like.0/sum(like.0, na.rm = TRUE)
        wts.0[is.na(wts.0)] <- 0
        return(list(wts.0 = wts.0, log.like.0 = log.like.0))
}



#' @title
#' Optimizer step for estimating the Heligman-Pollard Parameters using the Bayesian Melding with IMIS-opt procedure.
#'
#' @description
#' Performs the optimizer step in the IMIS procedure for the eight Heligman-Pollard parameters.
#'
#' @param prior A matrix containing the prior.
#' @param nrisk A vector containing the number of persons at risk in each age group.
#' @param ndeath A vector containing the number of deaths in each age group.
#' @param d Number of optimizer iterations.
#' @param theta.dim Number of columns of the prior (This should be 9 if estimating all parameters. Functionality for estimation a limited number of parameters does not exist yet).
#' @param age A vector containing the ages at which each age interval begins.
#'
#' @return
#' \item{opt.mu.d}{A matrix containing the local optimums resulting from the optimizer step. Each local optimum contains a set of 9 parameter values.}
#' \item{opt.cov.d}{An array containing the covariance matrix for each of the local optimums.}
#' \item{d.keep}{The number of local optimums found whose likelihood is greater than the maximum likelihood from the prior.}
#' \item{theta.new}{The set of parameters from the prior with the greatest weight as calculated with prior.likewts.}
#' \item{log.like.0}{A vector containing a likelihood for each row of the prior.}
#' \item{wts.0}{A vector containing an importance weight for each row of the prior.}
#'
#' @section Warning:
#' If the likelihood for the initial local maximum does not exceed the highlest likelihood from the prior, a warning will be issued.
#'
#' @note
#' Occasionally, this step fails to produce an initial local maximum that exceeds the highest likelihood of the prior and a warning is issued. Usually drawing a new prior or selecting a different algorithm solves this problem.
#'
#' @references
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{mod} \link{ll.binom} \link{prior.likewts} \link{HP.mod}
#'
#' @keywords Heligman-Pollard mortality bycatch
#'
#' @importFrom stats cov
#' @importFrom numDeriv hessian
#' @importFrom corpcor is.positive.definite
#' @importFrom numDeriv grad
#' @importFrom stats mahalanobis
#' @export

loop.optim <- function (prior, nrisk, ndeath, age, d = 10,
                        theta.dim = 9){
        lx <- nrisk
        dx <- ndeath
        H.k <- prior
        pllwts <- prior.likewts(prior = prior, nrisk = lx, ndeath = dx, age = age)
        log.like.0 <- pllwts$log.like.0
        wts.0 <- pllwts$wts.0
        B0 <- 1000 * theta.dim
        q0 <- H.k
        d.keep <- 0
        theta.new <- H.k[wts.0 == max(wts.0), ]
        keep <- H.k
        ll.keep <- log.like.0
        opt.mu.d <- matrix(NA, nrow = d, ncol = theta.dim)
        opt.cov.d <- array(NA, dim = c(theta.dim, theta.dim, d))
        prior.cov <- cov(q0)
        opt.low <- apply(q0, 2, min)
        opt.hi <- apply(q0, 2, max)
        imp.keep <- theta.dim * 100
        max.log.like.0 <- max(log.like.0)
        mp.mle <- function(theta, x.fit = age) {
                p.hat <- mod(theta = theta, x = x.fit)
                ll <- ll.binom(x = dx, n = lx, p = p.hat)
                return(ll)
        }
        for (i in 1:d) {
                out <- optim(par = theta.new, fn = mp.mle, method = "L-BFGS-B",
                             lower = opt.low, upper = opt.hi, control = list(fnscale = -1,
                                                                             maxit = 1e+05))
                out.mu <- out$par
                if (out$value > max.log.like.0) {
                        d.keep <- d.keep + 1
                        opt.mu.d[i, ] <- out.mu
                        out.hess <- hessian(func = mp.mle, x = out$par)
                        if (is.positive.definite(-out.hess)) {
                                out.cov <- try(solve(-out.hess))
                                opt.cov.d[, , i] <- out.cov
                        }
                        if (!is.positive.definite(-out.hess)) {
                                out.grad <- grad(func = mp.mle, x = out.mu)
                                A <- out.grad %*% t(out.grad)
                                out.prec <- try(solve(prior.cov)) + A
                                if (!is.positive.definite(out.prec)) {
                                        out.prec <- solve(prior.cov)
                                }
                                out.cov <- try(solve(out.prec))
                                opt.cov.d[, , i] <- out.cov
                        }
                }
                if (i == 1 & out$value <= max.log.like.0) {
                        out.hess <- hessian(func = mp.mle, x = out$par)
                        if (is.positive.definite(-out.hess)) {
                                out.cov <- solve(-out.hess)
                        }
                        if (!is.positive.definite(-out.hess)) {
                                out.grad <- grad(func = mp.mle, x = out.mu)
                                A <- out.grad %*% t(out.grad)
                                out.prec <- solve(prior.cov) + A
                                if (!is.positive.definite(out.prec)) {
                                        out.prec <- solve(prior.cov)
                                }
                                out.cov <- solve(out.prec)
                        }
                        warning("likelihood of first local maximum does not exceed maximum \t\t\tlikelihood from the prior")
                }
                if (i < d) {
                        keep <- keep[ll.keep != max(ll.keep), ]
                        ll.keep <- ll.keep[ll.keep != max(ll.keep)]
                        dist.to.mu <- mahalanobis(x = keep, center = out.mu,
                                                  cov = out.cov)
                        keep <- keep[rank(1/dist.to.mu) <= (d - i) * B0/d,
                                     ]
                        ll.keep <- ll.keep[rank(1/dist.to.mu) <= (d - i) *
                                                   B0/d]
                        theta.new <- keep[ll.keep == max(ll.keep), ]
                }
        }
        return(list(opt.mu.d = opt.mu.d, opt.cov.d = opt.cov.d, theta.new = theta.new,
                    d.keep = d.keep, log.like.0 = log.like.0, wts.0 = wts.0))
}



#' @title
#' Multivariate Gaussian Sampling for Heligman-Pollard model estimated via Bayesian Melding.
#'
#' @description
#' Samples the nine Heligman-Pollard parameters from the mvnorm distribution for each run of optimizer step where the likelihood for that run exceeds the maximum likelihood from the prior.
#'
#' @param opt.cov.d An array containing a covariance matrix for each run of optimizer where the likelihood for that run exceeds the maximum likelihood from the prior.
#' @param opt.mu.d A matrix containing the results of the optimizer step.
#' @param d.keep Number of runs of optimizer where the likelihood for that run exceeds the maximum likelihood from the prior.
#' @param prior A matrix containing the prior distribution (see \link{HP.priors}).
#' @param B Sample size at the importance sampling stage.
#' @param B0 Sample size of the prior. This is equal to (theta.dim * 1000).
#' @param d Number of optimizer iterations.
#'
#' @note
#' For use within the function \link{HP.mod}.
#'
#' @return
#' \item{H.k}{The prior plus new samples.}
#' \item{H.new}{The new samples from the multivariate normal.}
#' \item{B1}{The number of new samples - should be equal to B * d.keep.}
#'
#' @references
#'
#' Heligman, L. and Pollard, J.H. (1980) The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49-80.
#'
#' Poole, D and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244-1255.
#'
#' Raftery, A and Bao, L. (2009). "Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling." Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @keywords Heligman-Pollard post-sampling mortality bycatch
#'
#' @importFrom MASS mvrnorm
#' @export

samp.postopt <- function (opt.cov.d, opt.mu.d, d.keep, prior,
                          B = 400, B0 = 8000, d = 10)
{
        lo <- apply(prior, 2, min)
        hi <- apply(prior, 2, max)
        H.k <- prior
        for (i in 1:d) {
                if (!is.na(opt.mu.d[i, 1])) {
                        H.new <- mvrnorm(n = B, mu = opt.mu.d[i, ], Sigma = opt.cov.d[,
                                                                                      , i])
                        for (j in 1:B) {
                                lt0 <- table(H.new[j, ] <= 0)
                                while (as.numeric(lt0[1]) != ncol(prior) | H.new[j, 1] < lo[1] | H.new[j, 1] > hi[1] |
                                       H.new[j, 2] < lo[2] | H.new[j, 2] > hi[2] |
                                       H.new[j, 3] < lo[3] | H.new[j, 3] > hi[3] |
                                       H.new[j, 4] < lo[4] | H.new[j, 4] > hi[4] |
                                       H.new[j, 5] < lo[5] | H.new[j, 5] > hi[5] |
                                       H.new[j, 6] < lo[6] | H.new[j, 6] > hi[6] |
                                       H.new[j, 7] < lo[7] | H.new[j, 7] > hi[7] |
                                       H.new[j, 8] < lo[8] | H.new[j, 8] > hi[8] |
                                       H.new[j, 9] < lo[9] | H.new[j, 9] > hi[9]) {
                                        H.new[j, ] <- mvrnorm(1, mu = opt.mu.d[i, ],
                                                              Sigma = opt.cov.d[, , i])
                                        lt0 <- table(H.new[j, ] < 0)
                                }
                        }
                        H.k <- rbind(H.k, H.new)
                }
        }
        H.k <- cbind(H.k[, 1:9])
        H.new <- H.k[(B0 + 1):nrow(H.k), ]
        B1 <- d.keep * B
        return(list(H.k = H.k, H.new = H.new, B1 = B1))
}



#' @title
#' Local Optimums and Covariance from the optimizer step.
#'
#' @description
#' Defines some necessary arguments for the function \link{final.resamp}. Removes NAs from the opt.mu.d and opt.cov.d matrixes.
#'
#' @param K Number of iterations at the importance sampling stage.
#' @param log.like.0 A vector containing the likelihoods for each row of the prior.
#' @param opt.cov.d Covariance matrixes for the local optimums.
#' @param opt.mu.d A d x 8 matrix containing the local optimums (sets of parameters from the optimizer step).
#' @param d.keep Number of local optimums found in the optimizer step.
#' @param d A scalar defining the number of optimizer interations.
#' @param theta.dim Number of columns in the prior matrix.
#'
#' @return
#' \item{h.mu}{A d.keep * 8 matrix containing the local optimum result.}
#' \item{h.sig}{An array with (theta.dim * theta.dim * (K + d.keep)) dimensions containing the covariance matrix for each local optimum.}
#' \item{log.like}{A vector of likelihoods for each row of H.k.}
#'
#' @note
#' Typically for use immediately before running \link{final.resamp} or within the function \link{HP.mod}
#'
#' @references
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{final.resamp} \link{HP.mod}
#'
#' @keywords Heligman-Pollard likelihood resampling mortality bycatch
#'
#' @export

like.resamp <- function (K, log.like.0, opt.cov.d, opt.mu.d, d.keep, d = 10,
                         theta.dim = 9)
{
        K <- K
        log.like <- NULL
        log.like.k <- log.like.0
        h.mu <- NULL
        for (i in 1:d) {
                if (!is.na(opt.cov.d[1, 1, i])) {
                        h.mu <- rbind(h.mu, opt.mu.d[i, ])
                }
        }
        h.sig <- array(NA, dim = c(theta.dim, theta.dim, (K + d.keep)))
        m <- 0
        for (i in 1:d) {
                if (!is.na(opt.cov.d[1, 1, i])) {
                        m <- m + 1
                        h.sig[, , m] <- opt.cov.d[, , i]
                }
        }
        log.like <- log.like.k
        return(list(h.mu = h.mu, h.sig = h.sig, log.like = log.like,
                    K = K))
}



#' @title
#' Density of priors.
#'
#' @description
#' This function calculates the density of the prior distribution for the 9 parameters of the adapted Heligman-Pollard model. The density is calculated using a uniform distribution.
#'
#' @param x A 1 * 9 vector or n * 9 matrix containing values for the eight Heligman-Pollard Parameters.
#' @param pri.lo A vector giving the lower bounds of the uniform priors.
#' @param pri.hi A vector giving the upper bounds of the uniform priors.
#'
#' @return
#' A scalar describing the density of the prior distribution.
#'
#' @references
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @keywords Heligman-Pollard density priors mortality bycatch
#'
#' @importFrom stats dunif
#' @export

dens.prior <- function (x, pri.lo = c(0, 0, 0, 0, 0, 0, 15, 0, 0),
                        pri.hi = c(0.15, 1, 1, 0.5, 0.25, 15, 55, 0.1, 1.25))
{
        y <- (dunif(x[, 1], pri.lo[1], pri.hi[1]) *
                      dunif(x[, 2], pri.lo[2], pri.hi[2]) *
                      dunif(x[, 3], pri.lo[3], pri.hi[3]) *
                      dunif(x[, 4], pri.lo[4], pri.hi[4]) *
                      dunif(x[, 5], pri.lo[5], pri.hi[5]) *
                      dunif(x[, 6], pri.lo[6], pri.hi[6]) *
                      dunif(x[, 7], pri.lo[7], pri.hi[7]) *
                      dunif(x[, 8], pri.lo[8], pri.hi[8]) *
                      dunif(x[, 9], pri.lo[9], pri.hi[9]))
        return(y)
}



#' @title
#' Variance of the rescaled weights when estimating the Heligman-Pollard parameters using Bayesian Melding with IMIS.
#'
#' @description
#' Calculates the variance of the rescaled weights.
#'
#' @param w A vector of importance weights corresponding to each row of the mixture of the prior and multivariate gaussian draws.
#'
#' @return
#' A scalar representing the variance of the rescaled weights.
#'
#' @note
#' Used in the \link{final.resamp} function.
#'
#' @references
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{final.resamp} \link{HP.mod}
#'
#' @keywords Heligman-Pollard weights variance mortality bycatch
#'
#' @export

var.rwts <- function (w) {
        n <- length(w)
        ans <- mean((n * w - 1)^2)
        return(ans)
}



#' @title
#' Entropy of the rescaled weights relative to uniformity.
#'
#' @description
#' Performance measure for the IMIS algorithm that calculates the entropy of the importance weights relative to uniformity.
#'
#' @param w A vector of importance weights corresponding to each row of the mixture of the prior and multivariate gaussian draws.
#'
#' @return
#' Vector of entropy values relative to uniformity for a vector of weights.
#'
#' @note
#' For use in the function \link{final.resamp}.
#'
#' @references
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{final.resamp}
#'
#' @keywords Heligman-Pollard entropy mortality bycatch
#'
#' @export

entropy.wts <- function(w) {
        n <- length(w)
        num <- w * log(w)
        den <- log(n)
        ans <- -sum(num/den, na.rm = TRUE)
        return(ans)
}



#' @title
#' Expected number of unique inputs after the final IMIS re-sample.
#'
#' @description
#' Performance measure for the IMIS algorithm that calculates the expected number of unique points after re-sampling.
#'
#' @param w A vector of importance weights corresponding to each row of the mixture of the prior and multivariate gaussian distributions.
#' @param m The final re-sample size.
#'
#' @return
#' A scalar describing the number of unique points from the final re-sample.
#'
#' @note
#' For use in the function \link{final.resamp}.
#'
#' @references
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{final.resamp}
#'
#' @keywords Heligman-Pollard updates mortality bycatch
#'
#' @export

expt.upts <- function(w, m) {
        sum(1 - (1 - w)^m)
}



#' @title
#' Final re-sampling step in Bayesian Melding using IMIS.
#'
#' @description
#' Performs the final re-sampling step in the Bayesian Melding with IMIS procedure for the nine adapted Heligman-Pollard parameters.
#'
#' @param K The number of iterations of the importance sampling stage.
#' @param B1 Sample size at the importance sampling stage multiplied by the number of local optimums.
#' @param H.new A matrix with dimensions ((B * d.keep) * 9) containing the B*d.keep inputs drawn from the multivariate gaussians.
#' @param H.k A matrix containing the prior plus new inputs from the multivariate gaussians.
#' @param log.like A vector of log-likelihoods corresponding to each row of H.k.
#' @param d.keep The number of local optimums found in the optimizer step.
#' @param prior A matrix containing the prior.
#' @param h.mu A (d.keep * 9) matrix containing the results of the optimizer step.
#' @param h.sig An array containing the covariance matrix for each row of h.mu.
#' @param nrisk A vector containing the number of persons at risk in each age group.
#' @param ndeath A vector containing the number of deaths in each age group.
#' @param B Sample size at the importance sampling stage.
#' @param theta.dim The number of columns of the prior matrix.
#' @param age A vector containing the ages at which each age interval begins.
#'
#' @note
#' The function \link{HP.mod} performs this along with all other steps in a single function.
#' The algorithm ends when the expected fraction of unique points in the resample is at least 1 ??? 1/e = 0.632
#'
#' @return
#' \item{H.new}{A (B * theta.dim) matrix containing the posterior distribution for each parameter.}
#' \item{vwts}{A vector containing the variance of the rescaled weights at each IMIS iteration}
#' \item{ewts}{A vector containing the entropy of the rescaled weights at each IMIS iteration.}
#' \item{mwts}{A vector containing the maximum of the rescaled weights at each IMIS iteration.}
#' \item{nup}{A vector containing the expected number of unique points at each IMIS iteration.}
#' \item{frac.up}{A vector containing the proportion of unique points in the final resample at each IMIS iteration.}
#' \item{wts.k}{A vector containing the importance weights for the final iteration.}
#' \item{mwt.case}{The maximum weight value and associated case.}
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980) The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{mod} \link{ll.binom} \link{dens.prior} \link{var.rwts} \link{entropy.wts} \link{expt.upts}
#'
#' @keywords Heligman-Pollard resampling mortality bycatch
#'
#' @importFrom stats cov
#' @importFrom stats cov.wt
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats mahalanobis
#' @importFrom MASS mvrnorm
#' @export

final.resamp <- function (K, B1, H.new, H.k, log.like, d.keep, prior, h.mu, h.sig,
                          nrisk, ndeath, age, B = 400, theta.dim = 9) {
        lx <- nrisk
        dx <- ndeath
        mp.mle <- function(theta, x.fit = age) {
                p.hat <- mod(theta = theta, x = x.fit)
                ll <- ll.binom(x = dx, n = lx, p = p.hat)
                return(ll)
        }
        B0 <- theta.dim * 1000
        prior.cov <- cov(prior)
        imp.keep <- theta.dim * 100
        vwts <- NULL
        ewts <- NULL
        mwts <- NULL
        nup <- NULL
        nck <- NULL
        mwt.case <- NULL
        frac.nup <- NULL
        prilo <- apply(prior, 2, min)
        prihi <- apply(prior, 2, max)
        lo <- prilo
        hi <- prihi
        for (k in 1:K) {
                if (k == 1) {
                        log.like.k <- rep(NA, B1)
                        for (i in 1:B1) {
                                log.like.k[i] <- mp.mle(H.new[i, ])
                                if (i%%100 == 0) {
                                }
                        }
                }
                if (k > 1) {
                        log.like.k <- rep(NA, B)
                        for (i in 1:B) {
                                log.like.k[i] <- mp.mle(H.new[i, ])
                                if (i%%100 == 0) {
                                }
                        }
                        H.k <- rbind(H.k, H.new)
                }
                log.like <- c(log.like, log.like.k)
                n.k <- B * (k - 1) + B0 + B1
                ncol.qk <- 1 + d.keep + (k - 1)
                DENS.MAT <- matrix(NA, nrow = n.k, ncol = ncol.qk)
                DENS.MAT[, 1] <- (B0/n.k) * dens.prior(H.k, pri.lo = prilo, pri.hi = prihi)
                for (i in 1:(ncol.qk - 1)) {
                        DENS.MAT[, (i + 1)] <- B/n.k * dmvnorm(x = H.k, mean = h.mu[i,
                                                                                    ], sigma = h.sig[, , i])
                }
                qk <- rowSums(DENS.MAT)
                a <- -max(log.like, na.rm = TRUE)
                like.k <- exp(log.like)
                wts.k.0 <- like.k * dens.prior(H.k, pri.lo = prilo, pri.hi = prihi)
                wts.k <- wts.k.0/qk
                wts.k <- wts.k/sum(wts.k, na.rm = TRUE)
                wts.k[is.na(wts.k)] <- 0
                maxw.k <- max(wts.k)
                thresh.k <- min(maxw.k, 0.05)
                imp.theta <- H.k[wts.k == maxw.k, ]
                dist.to.mean <- mahalanobis(x = H.k, center = imp.theta,
                                            cov = prior.cov)
                keep.4.cov <- sort(dist.to.mean, decreasing = FALSE)[imp.keep]
                inputs.4.cov <- H.k[dist.to.mean <= keep.4.cov, ]
                wts.4.cov <- wts.k[dist.to.mean <= keep.4.cov]
                wts.4.cov <- ((wts.4.cov + 1/length(wts.k)) * 0.5)/sum(((wts.4.cov +
                                                                                 1/length(wts.k)) * 0.5))
                Sigma.k <- cov.wt(inputs.4.cov, wt = wts.4.cov)$cov
                h.mu <- rbind(h.mu, imp.theta)
                h.sig[, , (d.keep + k)] <- Sigma.k
                H.new <- mvrnorm(n = B, mu = imp.theta, Sigma = Sigma.k)
                for (i in 1:B) {
                        lt0 <- table(H.new[i, ] < 0)
                        while (lt0[1] != ncol(prior) |
                               H.new[i, 1] < lo[1] | H.new[i, 1] > hi[1] |
                               H.new[i, 2] < lo[2] | H.new[i, 2] > hi[2] |
                               H.new[i, 3] < lo[3] | H.new[i, 3] > hi[3] |
                               H.new[i, 4] < lo[4] | H.new[i, 4] > hi[4] |
                               H.new[i, 5] < lo[5] | H.new[i, 5] > hi[5] |
                               H.new[i, 6] < lo[6] | H.new[i, 6] > hi[6] |
                               H.new[i, 7] < lo[7] | H.new[i, 7] > hi[7] |
                               H.new[i, 8] < lo[8] | H.new[i, 8] > hi[8] |
                               H.new[i, 9] < lo[9] | H.new[i, 9] > hi[9]) {
                                H.new[i, ] <- mvrnorm(n = 1, mu = imp.theta,
                                                      Sigma = Sigma.k)
                                lt0 <- table(H.new[i, ] < 0)
                        }
                        if (i%%100 == 0) {
                        }
                }
                H.new <- cbind(H.new[, 1:9])
                vwts <- c(vwts, var.rwts(wts.k))
                ewts <- c(ewts, entropy.wts(wts.k))
                mwts <- c(mwts, maxw.k)
                nup <- c(nup, expt.upts(wts.k, m = B))
                frac.nup <- c(frac.nup, expt.upts(wts.k, m = B)/B)
                N.k <- K * B
                case.n <- 1:N.k
                mwt.case <- c(mwt.case, case.n[wts.k == maxw.k])
                assign(paste("imp.theta.k", k, sep = ""), imp.theta)
                assign(paste("Sigma.k", k, sep = ""), Sigma.k)
                if (k%%5 == 0) {
                        assign(paste("log.like.k", k, sep = ""), log.like.k)
                        assign(paste("wts.k", k, sep = ""), wts.k)
                        assign(paste("maxw.k", k, sep = ""), maxw.k)
                        assign(paste("imp.theta", k, sep = ""), imp.theta)
                        assign(paste("H.k", k, sep = ""), H.k)
                        assign(paste("qk", k, sep = ""), qk)
                        save.vectors <- c(paste("log.like.k", k, sep = ""),
                                          paste("wts.k", k, sep = ""), paste("maxw.k",
                                                                             k, sep = ""), paste("imp.theta", k, sep = ""),
                                          paste("Sigma.k", k, sep = ""), paste("H.k", k,
                                                                               sep = ""), "qk", "log.like")
                        rm.vectors <- c(paste("log.like.k", k, sep = ""),
                                        paste("wts.k", k, sep = ""), paste("maxw.k",
                                                                           k, sep = ""))
                        rm(list = rm.vectors)
                }
                if (frac.nup[k] >= 0.632) {
                        break
                }
        }
        return(list(H.new = H.new, vwts = vwts, ewts = ewts, mwts = mwts,
                    nup = nup, frac.nup = frac.nup, wts.k = wts.k, mwt.case = mwt.case))
}



#' @title
#' Heligman-Pollard parameter estimator using Bayesian Melding with Incremental Mixture Importance Sampling.
#'
#' @description
#' Runs all the necessary functions to estimate the nine adapted Heligman-Pollard parameters in one step via
#' Bayesian Melding with IMIS and optimization. In this order and with the proper arguments imputed
#' the functions run are \link{loop.optim}, \link{samp.postopt}, \link{like.resamp}, \link{final.resamp}.
#'
#' @param prior A matrix with dimensions (9000 * theta.dim) containing the prior distribution for each Heligman-Pollard parameter.
#' @param lifeTab A life table object from \link{life.tab} with a column containing the total number of individuals at risk of death in each age group.
#' A column containing the total number of deaths in each age group. Length should equal the length of age.
#' and a column of the ages at which the probabilities of death will be calculated.
#' @param K The number of IMIS iterations.
#' @param d The number of optimizer iterations.
#' @param B The sample size at each importance sampling iteration.
#' @param CI Defines the width of the confidence interval.
#' A summary table is printed with the median estimate and lower and upper confidence bounds.
#' Setting CI=95 prints a table with the first column representing the 2.5th percentile for each parameter distribution,
#' the second column represents the median value for each parameter distribution
#' and the third column represents the 97.5th percentile for each parameter distribution.
#' @param rm The number of age classes that want to be removed from optimization
#'
#' @details
#' All age classes are used for adjustment, in case of whish to remove any of the first age classes due to bias in the sample, indicate it with the parameter "rm" and these will be removed starting from the first age class
#'
#' @return
#' A list with:
#' \item{out}{A summary table of the results with the median parameter values in the middle column, the lower bound results in the left column, and upper bound result in the right column.}
#' \item{H.final}{A (B * theta.dim) matrix containing the posterior distribution for each parameter.}
#' \item{h.mu}{The sets of parameters found in the optimizer step.}
#' \item{h.sig}{The covariance matrix for each set of parameters in h.mu.}
#' \item{log.like}{A vector of likelihoods for the prior plus resamples.}
#' \item{log.like.0}{A vector of the likelihoods for the prior.}
#' \item{wts.0}{A vector of importance weights for each set of parameters in the prior.}
#' \item{d.keep}{The number of optimizer runs where the likelihood exceeded the maximum likelihood of the prior.}
#' \item{vwts}{A vector containing the variance of the rescaled weights at each IMIS iteration.}
#' \item{ewts}{A vector containing the entropy of the rescaled weights at each IMIS iteration.}
#' \item{mwts}{A vector containing the maximum of the rescaled weights at each IMIS iteration.}
#' \item{mwt.case}{The maximum weight and associated case.}
#' \item{nup}{A vector containing the expected number of unique points at each IMIS iteration.}
#' \item{frac.up}{A vector containing the proportion of unique points in the final resample at each IMIS iteration.}
#' \item{wts.k}{A vector containing the importance weights for the final IMIS iteration.}
#'
#' @note
#' Because there are multiple sampling steps sometimes with upper and lower bound restricitions, this function can take several minutes to run depending on the sample size, K.
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980) The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{loop.optim} \link{samp.postopt} \link{like.resamp} \link{final.resamp}
#'
#' @keywords Heligman-Pollard mortality bycatch
#'
#' @importFrom stats quantile
#' @importFrom stats median
#' @export

HP.mod <- function (prior, lifeTab, K = 10, d = 10, B = 500, CI = 95, rm = 0) {
        if(rm == 0){
                nrisk = round(lifeTab$lx * 1000)
                ndeath = round(lifeTab$dx)
                age = lifeTab$age
        } else {
                nrisk = round(lifeTab$lx * 1000)[-c(1:rm)]
                ndeath = round(lifeTab$dx)[-c(1:rm)]
                age = lifeTab$age[-c(1:rm)]
        }
        low <- apply(prior, 2, min)
        high <- apply(prior, 2, min)
        opt.result <- loop.optim(prior = prior, nrisk = nrisk, ndeath = ndeath,
                                 d = d, theta.dim = ncol(prior), age = age)
        opt.mu.d <- opt.result$opt.mu.d
        opt.cov.d <- opt.result$opt.cov.d
        theta.new <- opt.result$theta.new
        d.keep <- opt.result$d.keep
        log.like.0 <- opt.result$log.like.0
        wts.0 <- opt.result$wts.0
        samp.po <- samp.postopt(opt.cov.d = opt.cov.d, opt.mu.d = opt.mu.d,
                                prior = prior, d.keep = d.keep, B = B, d = d, B0 = nrow(prior))
        H.k <- samp.po$H.k
        H.new <- samp.po$H.new
        B1 <- samp.po$B1
        ll.postopt <- like.resamp(K = K, log.like.0 = log.like.0,
                                  opt.cov.d = opt.cov.d, opt.mu.d = opt.mu.d, d.keep = d.keep,
                                  d = d, theta.dim = ncol(prior))
        h.mu <- ll.postopt$h.mu
        h.sig <- ll.postopt$h.sig
        log.like <- ll.postopt$log.like
        K <- ll.postopt$K
        result <- final.resamp(K = K, B1 = B1, B = B, H.new = H.new,
                               H.k = H.k, log.like = log.like, d.keep = d.keep, prior = prior,
                               h.mu = h.mu, h.sig = h.sig, nrisk = nrisk, ndeath = ndeath,
                               theta.dim = ncol(prior), age = age)
        H.final <- result$H.new
        nup <- result$nup
        vwts <- result$vwts
        ewts <- result$ewts
        frac.nup <- result$frac.nup
        mwts <- result$mwts
        wts.k <- result$wts.k
        mwt.case <- result$mwt.case
        loCI <- ((100 - CI)/2)/100
        hiCI <- 1 - (((100 - CI)/2)/100)
        loci <- round(rbind(quantile(H.final[, 1], loCI),
                            quantile(H.final[, 2], loCI),
                            quantile(H.final[, 3], loCI),
                            quantile(H.final[, 4], loCI),
                            quantile(H.final[, 5], loCI),
                            quantile(H.final[, 6], loCI),
                            quantile(H.final[, 7], loCI),
                            quantile(H.final[, 8], loCI),
                            quantile(H.final[, 9], loCI)), digits = 3)
        Median <- round(rbind(median(H.final[, 1]),
                              median(H.final[, 2]),
                              median(H.final[, 3]),
                              median(H.final[, 4]),
                              median(H.final[, 5]),
                              median(H.final[, 6]),
                              median(H.final[, 7]),
                              median(H.final[, 8]),
                              median(H.final[, 9])), digits = 3)
        hici <- round(rbind(quantile(H.final[, 1], hiCI),
                            quantile(H.final[, 2], hiCI),
                            quantile(H.final[, 3], hiCI),
                            quantile(H.final[, 4], hiCI),
                            quantile(H.final[, 5], hiCI),
                            quantile(H.final[, 6], hiCI),
                            quantile(H.final[, 7], hiCI),
                            quantile(H.final[, 8], hiCI),
                            quantile(H.final[, 9], hiCI)), digits = 3)
        pout <- data.frame(loci, Median, hici)
        names(pout) <- c("Low CI", "Median", "High CI")
        print(pout)
        return(list(H.final = H.final, h.mu = h.mu, h.sig = h.sig,
                    log.like = log.like, log.like.0 = log.like.0, wts.0 = wts.0,
                    d.keep = d.keep, nup = nup, frac.nup = frac.nup, vwts = vwts,
                    ewts = ewts, mwts = mwts, wts.k = wts.k, mwt.case = mwt.case,
                    out = pout))
}

