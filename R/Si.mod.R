#' @title
#' Siler model.
#'
#' @description
#' Fit a 5-parameters Competing-Risk Siler model for Animal Mortality.
#'
#' @param data Data frame with age clases and frequency of occurrence (see Details).
#' @param par Initial values for the Siler parameters to be optimized over.
#' @param rm The number of age classes that want to be removed from optimization (see Details).
#' @param method The method to be used: "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"
#'   or "Brent" (see \link{optim}).
#' @param control A list of control parameters (see \link{optim}).
#'
#' @details
#' The data used must be a data frame whose first column is a vector of the estimated ages of the animals found dead and the second the frequency of occurrence of those ages.
#'
#' All age classes are used for adjustment, in case of whish to remove any of the first age classes due to bias in the sample, indicate it with the parameter "rm" and these will be removed starting from the first age class.
#'
#' @references
#' Siler, W. (1979). A Competing-Risk Model for Animal Mortality. Ecology 60, 750–757.
#'
#' Siler, W. (1983). Parameters of mortality in human populations with widely varying life spans. Stat. Med. 2, 373–380.
#'
#' Nocedal, J. and Wright, S. J. (1999). Numerical Optimization. Springer.
#'
#' @seealso \link{optim}
#'
#' @keywords Siler mortality
#'
#' @importFrom stats optim
#'
#' @examples
#'
#' Si.mod(data = cetaceans, rm = 2,
#'        par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))
#'
#' Si.mod(data = cetaceans, rm = 1,
#'        par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))
#'
#' Si.mod(data = cetaceans, rm = 0,
#'        par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))
#'
#' @export


Si.mod <- function(data, par = c(-0.15, 1.10, 0.15, 0.005, 0.15),
                   rm = 0, method = "Nelder-Mead",
                   control = list(fnscale = -1, maxit = 10000)){
        data$age1 <- data[,1] + 1
        optim(par,
              function(par) {
                      a1 <- par[1]
                      b1 <- par[2]
                      a2 <- par[3]
                      a3 <- par[4]
                      b3 <- par[5]
                      if(rm != 0){
                              data <- data[-c(1,rm),]
                              data[,1] <- data[,1] - rm
                              data[,length(data)] <- data[,length(data)] - rm
                      }
                      n <- NROW(data)
                      S.t <- function(t) {
                              return(exp(-a1/b1 * (1 - exp(-b1 * t))) *
                                             exp(-a2 * t) * exp(a3/b3 * (1 - exp(b3 * t))))
                      }
                      dif <- S.t(data[1:n,1]) - S.t(data[1:n,length(data)])
                      obs <- data[,2]
                      lnlk <- as.numeric(crossprod(obs, log(dif)))
                      return(lnlk)
              },
              method = method,
              control = control)
}
