#' @title
#' Estimation of starting values for priors of the Heligman-Pollard model.
#'
#' @description
#' Estimates starting values for the priors of the modified 9-parameter Heligman-Pollard model from mortality rates by age and proportion of individuals death from non natural causes.
#'
#' @param Si.qx A vector of the same length as "age" containing the total mortality rates by age (qx).
#' @param Prop.risk A vector of the same length as "age" containing the proportion of death animals identified as death by external effects (other different than natural mortality).
#' @param Life.qx Mortality vector qx form life table created with \link{life.tab}.
#' @param age A vector containing the ages at which each age interval begins.
#' @param rg The variation of the range at which the mean prior values will be multiplied to calculate low and high limits.
#' @param control Allow the user to set some characteristics Levenberg-Marquardt nonlinear least squares algorithm implemented in \link{nls.lm}. see \link{nls.lm.control}
#'
#' @return
#' A dataframe with low, high and mean starting values for priors, to be used for fitting a adapted 9-parametres Heligman-Pollard model
#'
#' @references
#' Heligman, L. and Pollard, J.H. (1980). The Age Pattern of Mortality. Journal of the Institute of Actuaries 107:49–80.
#'
#' Moré, J.J. (1978). The Levenberg-Marquardt algorithm: implementation and theory, in: Lecture Notes in Mathematics 630: Numerical Analysis, G.A. Watson (Ed.), Springer-Verlag: Berlin, pp. 105-116.
#'
#' @seealso \link{HP.priors} \link{HP.mod} \link{nlsLM} \link{life.tab}
#'
#' @keywords Heligman-Pollard priors mortality bycatch
#'
#' @importFrom minpack.lm nls.lm.control
#' @importFrom minpack.lm nlsLM
#' @importFrom stats coef
#' @export

HP.pri.start <- function(Si.qx, Prop.risk, Life.qx, age, rg = 0.25,
                            control = nls.lm.control(maxiter = 100)){
        Y <- Si.qx * (1 - Prop.risk)
        n <- which(Y[-length(Y)] - Y[-1] < 0)[1]
        Y[n:length(age)]<- Y[n]
        my <- nlsLM(Y ~ a^(((age) + b)^c),
                    start = list(a = 0.05, b = 0.5, c = 1),
                    control= control)
        R <- Life.qx * Prop.risk
        mr <- nlsLM(R ~ a + b * exp(-c * (log(age) - log(d))^2),
                    start = list(a = 0.05, b = 0.01, c = 0.1, d = 1),
                    control= control)
        A <- Si.qx * (1 - Prop.risk)
        nA <- which(A[-length(A)] - A[-1] < 0)[1]-1
        A[1:nA]<- Y[nA]
        ma <- nlsLM(A ~ (a * (b^(age))) / (1 + a * (b^(age))),
                    start =list( a = 0.001, b = 1.1),
                    control = control)
        start.priors <-
                data.frame(
                        priors.lo =
                                c(ifelse(coef(my)[1]<0, coef(my)[1]+coef(my)[1]*rg, coef(my)[1]-coef(my)[1]*rg),
                                  ifelse(coef(my)[2]<0, coef(my)[2]+coef(my)[2]*rg, coef(my)[2]-coef(my)[2]*rg),
                                  ifelse(coef(my)[3]<0, coef(my)[3]+coef(my)[3]*rg, coef(my)[3]-coef(my)[3]*rg),
                                  ifelse(coef(mr)[1]<0, coef(mr)[1]+coef(mr)[1]*rg, coef(mr)[1]-coef(mr)[1]*rg),
                                  ifelse(coef(mr)[2]<0, coef(mr)[2]+coef(mr)[2]*rg, coef(mr)[2]-coef(mr)[2]*rg),
                                  ifelse(coef(mr)[3]<0, coef(mr)[3]+coef(mr)[3]*rg, coef(mr)[3]-coef(mr)[3]*rg),
                                  ifelse(coef(mr)[4]<0, coef(mr)[4]+coef(mr)[4]*rg, coef(mr)[4]-coef(mr)[4]*rg),
                                  ifelse(coef(ma)[1]<0, coef(ma)[1]+coef(ma)[1]*rg, coef(ma)[1]-coef(ma)[1]*rg),
                                  ifelse(coef(ma)[2]<0, coef(ma)[2]+coef(ma)[2]*rg, coef(ma)[2]-coef(ma)[2]*rg)),
                        priors.med =
                                c(coef(my)[1], coef(my)[2], coef(my)[3],
                                  coef(mr)[1], coef(mr)[2], coef(mr)[3], coef(mr)[4],
                                  coef(ma)[1], coef(ma)[2]),
                        priors.hi =
                                c(ifelse(coef(my)[1]<0, coef(my)[1]-coef(my)[1]*rg, coef(my)[1]+coef(my)[1]*rg),
                                  ifelse(coef(my)[2]<0, coef(my)[2]-coef(my)[2]*rg, coef(my)[2]+coef(my)[2]*rg),
                                  ifelse(coef(my)[3]<0, coef(my)[3]-coef(my)[3]*rg, coef(my)[3]+coef(my)[3]*rg),
                                  ifelse(coef(mr)[1]<0, coef(mr)[1]-coef(mr)[1]*rg, coef(mr)[1]+coef(mr)[1]*rg),
                                  ifelse(coef(mr)[2]<0, coef(mr)[2]-coef(mr)[2]*rg, coef(mr)[2]+coef(mr)[2]*rg),
                                  ifelse(coef(mr)[3]<0, coef(mr)[3]-coef(mr)[3]*rg, coef(mr)[3]+coef(mr)[3]*rg),
                                  ifelse(coef(mr)[4]<0, coef(mr)[4]-coef(mr)[4]*rg, coef(mr)[4]+coef(mr)[4]*rg),
                                  ifelse(coef(ma)[1]<0, coef(ma)[1]-coef(ma)[1]*rg, coef(ma)[1]+coef(ma)[1]*rg),
                                  ifelse(coef(ma)[2]<0, coef(ma)[2]-coef(ma)[2]*rg, coef(ma)[2]+coef(ma)[2]*rg))
                )
        return(start.priors)
}
