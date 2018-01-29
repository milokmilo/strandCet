#' @title
#' Prediction of Heligman-Pollard model.
#'
#' @description
#' Predicts Heligman-Pollard model from Heligman-Pollard's parameters.
#'
#' @param life A life table created with \link{life.tab} or a dataframe with a vector of ages in the first column.
#' @param HPout A model object created with \link{HP.mod} with the Heligman-Pollard estimated params.
#' @param M Defines the statistic to predict. Median by default (med), Low CI (low) or High CI (high)
#' @param age A vector containing the ages at which each age interval begins. See Details
#' @param rm The number of age classes that want to be removed from optimization.
#'
#' @details
#' Mx is returned only if number of ages required for prediction is equal to the number of ages in the life table.
#'
#' @return
#' A dataframe with seven columns:
#' \item{age}{Age at the beginning of the interval.}
#' \item{Mx}{Number of observed deaths at age x.}
#' \item{qx.tot}{Total probability of death between ages x and x + n.}
#' \item{qx.nat}{Natural probability of death between ages x and x + n.}
#' \item{qx.young}{Young probability of death between ages x and x + n.}
#' \item{qx.risk}{Probability of death due to an externl risk between ages x and x + n.}
#' \item{qx.adult}{Adult or senescent probability of death between ages x and x + n.}
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
#'
#' @seealso \link{HP.mod}
#'
#' @keywords Heligman-Pollard mortality prediction
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
#' modHP <- HP.mod(prior = q0, lifeTab = lifeN, rm = 2, K = 10, d = 10, B = 10, CI = 90)
#'
#' HP.pred(life = lifeN, HPout = modHP, age = seq(0,29,1), rm = 2)
#'
#' @export

HP.pred <- function(life, HPout, M = "med", age = seq(0,29,0.1), rm = 0){
        if(M == "med"){
                A <- HPout$out[1, 2]
                B <- HPout$out[2, 2]
                C <- HPout$out[3, 2]
                D2 <- HPout$out[4, 2]
                D <- HPout$out[5, 2]
                E <- HPout$out[6, 2]
                F <- HPout$out[7, 2]
                G <- HPout$out[8, 2]
                H <- HPout$out[9, 2]
                }
        if(M == "high"){
                A <- HPout$out[1, 3]
                B <- HPout$out[2, 3]
                C <- HPout$out[3, 3]
                D2 <- HPout$out[4, 3]
                D <- HPout$out[5, 3]
                E <- HPout$out[6, 3]
                F <- HPout$out[7, 3]
                G <- HPout$out[8, 3]
                H <- HPout$out[9, 3]
                }
        if(M == "low"){
                A <- HPout$out[1, 1]
                B <- HPout$out[2, 1]
                C <- HPout$out[3, 1]
                D2 <- HPout$out[4, 1]
                D <- HPout$out[5, 1]
                E <- HPout$out[6, 1]
                F <- HPout$out[7, 1]
                G <- HPout$out[8, 1]
                H <- HPout$out[9, 1]
                }
        age <- age
        Myoung <- A^(((age) + B)^C)
        Mrisk <- D2 + D * exp(-E * (log(age) - log(F))^2)
        Madult <- (G * (H^(age))) / (1 + G * (H^(age)))
        Mtot <- Myoung + Mrisk + Madult
        Mnat <- Myoung + Madult
        if (length(age) == length(life$Mx)){
                Mx <- life$Mx
                if(!rm == 0){
                        Mx[1:rm] <- NA
                }
                dataHP <- data.frame(age = age, Mx = Mx, qx.tot = Mtot, qx.nat = Mnat,
                                     qx.young = Myoung, qx.risk = Mrisk, qx.adult = Madult)
        } else {
                dataHP <- data.frame(age = age, qx.tot = Mtot, qx.nat = Mnat,
                                     qx.young = Myoung, qx.risk = Mrisk, qx.adult = Madult)
        }
        return(dataHP)
}

