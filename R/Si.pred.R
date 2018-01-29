#' @title
#' Predict Siler model
#'
#' @description
#' Predict Siler model from Siler's parameters
#'
#' @param data frame with age clases and frequency of occurrence (see Details).
#' @param Sout A model object created with \link{Si.mod} with the 5 Siler's parameters.
#' @param rm The number of age classes that want to be removed from optimization (see Details).
#'
#' @details
#' The data used must be the data frame from which the Siler's patameters were estimated, whose first column is a vector of the estimated ages of the animals found dead and the second the frequency of occurrence of those ages.
#'
#' In case that any age class had been removed for the estimation of the parameters that are going to be used for predicting the model, indicate it with the parameter "rm" and these age classes will be removed starting from the first age class.
#'
#' @references
#' Siler, W. (1979). A Competing-Risk Model for Animal Mortality. Ecology 60, 750–757.
#'
#' Siler, W. (1983). Parameters of mortality in human populations with widely varying life spans. Stat. Med. 2, 373–380.
#'
#' @seealso \link{Si.mod}
#'
#' @keywords Siler mortality
#'
#' @examples
#'
#' modSi <- Si.mod(data = cetaceans, rm = 2,
#'                 par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))
#'
#' Si.pred(data = cetaceans, Sout = modSi, rm = 2)
#'
#' @export

Si.pred <- function(data, Sout, rm=0){
        a1 <- Sout$par[1]
        b1 <- Sout$par[2]
        a2 <- Sout$par[3]
        a3 <- Sout$par[4]
        b3 <- Sout$par[5]
        t <- data[,c(1,2)]
        if(!rm == 0){
                t[1:rm, 2] <- NA
        }
        age <- data$age - rm
        Mtot <- a1 * exp(-b1 * age) + a2 + a3 * exp(b3 * age)
        Mcte <- rep(a2, length(age))
        Myoung <- a1 * exp(-b1 * age)
        Madult <- a3 * exp(b3 * age)
        dataSiler <-  data.frame(t, qx.tot = Mtot,
                                 qx.young = Myoung, qx.cte = Mcte, qx.adult = Madult)
        return(dataSiler)
}

