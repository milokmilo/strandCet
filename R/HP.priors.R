#' @title
#' Heligman-Pollard Parameter prior formation.
#'
#' @description
#' Draws from a uniform distribution with bounds "pri.lo" and "pri.hi" to create the prior distribution
#' of the Heligman-Pollard parameters necessary for the Bayesian Melding procedure.
#'
#' @param pri.lo Lower bound of the uniform from which the prior is drawn.
#' @param pri.hi Upper bound of the uniform from which the prior is drawn.
#' @param theta.dim The number of parameters to be estimated.
#'
#' @return
#' A ((1000 * theta.dim) x theta.dim) matrix containing the 1000 * theta.dim sets of the Heligman-Pollard
#' parameters drawn from a uniform distribution.
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.
#'
#' Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @keywords Heligman-Pollard priors mortality bycatch
#'
#' @examples
#'
#' priors <- data.frame(priors.lo = c(0,0.5,0,0,0,0,6,0,1),
#'                      priors.hi = c(0.1,1,1,0.15,0.15,50,10,0.01,1.5))
#'
#' HP.priors(pri.lo = priors$priors.lo,
#'           pri.hi = priors$priors.hi,
#'           theta.dim = 9)
#'
#' @importFrom stats runif
#' @export

HP.priors <- function (pri.lo = c(0, 0, 0, 0.001, 0, 0, 15, 0, 0),
                       pri.hi = c(0.15, 1, 1, 0.5, 0.25, 15, 55, 0.1, 1.25),
                       theta.dim = 9) {
        B0 <- 1000 * theta.dim
        q0 <- cbind(runif(B0, pri.lo[1], pri.hi[1]),
                    runif(B0, pri.lo[2], pri.hi[2]),
                    runif(B0, pri.lo[3], pri.hi[3]),
                    runif(B0, pri.lo[4], pri.hi[4]),
                    runif(B0, pri.lo[5], pri.hi[5]),
                    runif(B0, pri.lo[6], pri.hi[6]),
                    runif(B0, pri.lo[7], pri.hi[7]),
                    runif(B0, pri.lo[8], pri.hi[8]),
                    runif(B0, pri.lo[9], pri.hi[9]))
        H.k <- q0
        return(H.k)
}

