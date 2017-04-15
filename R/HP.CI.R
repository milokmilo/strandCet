#' @title
#' Heligman-Pollard parameter coversion to age-specific probabilites of death due to an external risk.
#'
#' @description
#' Calculates the age-specific probabilities of death due to an external risk using the Heligman-Pollard model.
#'
#' @param theta A vector containing values for the 9 parameters of the adapted Heligman-Pollard model.
#' @param x A vector containing the ages at which to calculate the probabilities of death.
#'
#' @note
#' Utility function called by \code{\link{HP.CI}}.
#'
#' @return
#' A vector of probabilities of death at ages due to an external risk.
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' @seealso \code{\link{HP.CI}}
#'
#' @keywords Heligman-Pollard bycatch mortality
#'
#' @export

mod.risk <- function (theta, x) {
        D2 <- theta[4]
        D <- theta[5]
        E <- theta[6]
        F <- theta[7]
        f.x <- (D2 + D * exp(-E * (log(x) - log(F))^2))
        f.x2 <- f.x
        f.x2[1e-06 > f.x] <- 1e-06
        f.x2[0.999999 < f.x] <- 0.999999
        return(f.x2)
}



#' @title
#' Heligman-Pollard parameter coversion to natural age-specific probabilites of death.
#'
#' @description
#' Calculates the age-specific probabilities of death using the Heligman-Pollard model.
#'
#' @param theta A vector containing values for the 9 parameters of the adapted Heligman-Pollard model.
#' @param x A vector containing the ages at which to calculate the probabilities of death.
#'
#' @note
#' Utility function called by \code{\link{HP.CI}}.
#'
#' @return
#' A vector of probabilities of natural death at ages.
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' @seealso \code{\link{HP.CI}}
#'
#' @keywords Heligman-Pollard mortality bycatch
#'
#' @export

mod.nat <- function (theta, x) {
        A <- theta[1]
        B <- theta[2]
        C <- theta[3]
        G <- theta[8]
        H <- theta[9]
        f.x <- A^((x + B)^C) + (G * (H^x))/(1 + G * (H^x))
        f.x2 <- f.x
        f.x2[1e-06 > f.x] <- 1e-06
        f.x2[0.999999 < f.x] <- 0.999999
        return(f.x2)
}



#' @title
#' Heligman-Pollard parameter conversion to age-specifc probabilites of death.
#'
#' @description
#' Converts a set of Heligman-Pollard mortality model parameters into age-specific probabilities of death.
#'
#' @param HPout A model object created with \link{HP.mod} with the Heligman-Pollard estimated params.
#' @param age A vector containing the ages at which the probability of death will be calculated.
#' @param M The type of probabilities to be calculated. Can be both "natural" mortality, "total" mortality and mortality due to an external risk.
#'
#' @note
#' Utility function called by \code{\link{HP.CI}}.
#'
#' @return
#' Set of age specific probabilities of death equal to the length of age.
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' @seealso \code{\link{HP.CI}}
#'
#' @keywords Heligman-Pollard mortality probability
#'
#' @export

hp.nqx <- function (HPout, age, M = "total") {
        H.new.hat <- matrix(NA, nrow = nrow(HPout$H.final), ncol = length(age))
        for (i in 1:nrow(HPout$H.final)) {
                if(M == "total"){H.new.hat[i, ] <- mod(theta = HPout$H.final[i, ], x = age)}else{
                        if(M == "natural"){H.new.hat[i, ] <- mod.nat(theta = HPout$H.final[i, ], x = age)}else{
                                if(M == "risk"){H.new.hat[i, ] <- mod.risk(theta = HPout$H.final[i, ], x = age)}else{
                                        stop("M parameter wrong defined")}}}
        }
        ans <- H.new.hat
        return(ans)
}



#' @title
#' Helligman-Pollard confidence intervals with 9 parameters
#'
#' @description
#' Predicts Helligman-Pollard model with the estimated parameters both "natural" mortality,
#' "total" mortality and mortality due to an external risk with confidence intervals.
#'
#' @param HPout A model object created with \link{HP.mod} with the Heligman-Pollard estimated params
#' @param age A vector of the ages at which the probabilities of death will be calculated.
#' @param CI Defines the width of the credible interval (Defaults to 95 percent).
#' A summary table is printed with the median estimate and lower and upper confidence bounds.
#' Setting CI = 95 prints a table with the first column representing the 2.5th
#' percentile for each parameter distribution, the second column represents the median
#' value for each parameter distribution and the third column represents the
#' 97.5th percentile for each parameter distribution.
#' @param M The type of probabilities to be calculated. Can be both "natural" mortality,
#' "total" mortality and mortality due to an external risk.
#'
#' @details
#' The type of mortality to be calculated must be defined. By default the total probability of death
#' is calculated but only natural mortality or due to an external risk can be calculad if previously defined.
#'
#' @return
#' Return a dataframe with number of rows equal to the number of age classes and four columns:
#' \item{age}{Vetor with age classes.}
#' \item{Med}{Median prediction of probabilities of death at age.}
#' \item{Mlo}{Lower limit of prediction of probabilities of death at age with CI = (100 - CI)/2.}
#' \item{Mhi}{Higher limit of prediction of probabilities of death at age with CI = 1- (100 - CI)/2.}
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' @seealso \code{\link{hp.nqx}} \code{\link{mod}} \code{\link{mod.nat}} \code{\link{mod.risk}}
#'
#' @keywords Heligman-Pollard confidence-intervals
#'
#' @export


HP.CI <- function(HPout, age, CI = 95, M = "total"){
        loCI <- ((100 - CI)/2)/100
        hiCI <- 1 - (((100 - CI)/2)/100)
        hpq <- hp.nqx(HPout = HPout, age = age, M = M)
        hpq.med <- rep(NA, length(age))
        for (i in 1:length(hpq.med)) {
                hpq.med[i] <- median(hpq[, i])
        }
        hpq5 <- rep(NA, length(age))
        for (i in 1:length(hpq5)) {
                hpq5[i] <- quantile(hpq[, i], probs = loCI)
        }
        hpq95 <- rep(NA, length(age))
        for (i in 1:length(hpq95)) {
                hpq95[i] <- quantile(hpq[, i], probs = hiCI)
        }
        MedLoHi <- data.frame(age = age, Med = hpq.med, Mlo = hpq5, Mhi = hpq95)
        return(MedLoHi)
}
