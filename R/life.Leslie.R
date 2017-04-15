#' @title
#' Cohort method
#'
#' @description
#' Dummy function called in \code{\link{life.Leslie}}
#'
#' @param width12 Width of the first two age classes.
#'
#' @note
#' Utility function called by \code{\link{life.Leslie}}.
#'
#' @return
#' A two element vector representing the first two values of the nax column of the life table.
#'
#' @seealso \code{\link{life.Leslie}}
#'
#' @keywords cohort life-table leslie-matrix
#'
#' @export

cohort <- function(width12){
        return(nax12 <- width12/2)
}


#' @title
#' Coale method.
#'
#' @description
#' Used to graduate the individuals-years lived by those dying in the interval by the method of Coale et al. (1983).
#'
#' @param b1 Two-element vector of regression coefficients for graduating 1 to 0 age-classes provided in Coale et al. (1983).
#' @param b4 Two-element vector of regression coefficients for graduating 4 to 1 age-classes provided in Coale et al. (1983).
#' @param nMx Period central death rates = nDx/nKx.
#'
#' @note
#' Utility function called by \code{\link{life.Leslie}}.
#'
#' @keywords coale life-table leslie-matrix
#'
#' @export

coale <- function(b1,b4,nMx){
        if(nMx[1]>0.107){
                b1 <- c(0.350,0)
                b4 <- c(1.361,0)
        }
        nax12 <- c(0,0)
        nax12[1] <- b1[1] + b1[2] *nMx[1]
        nax12[2] <- b4[1] + b4[2]* nMx[1]
        return(nax12)
}


#' @title
#' Keyfitz and Flieger method.
#'
#' @description
#' Utility used by \code{\link{life.Leslie}} to graduate the person-years lived by those dying in the interval by the
#' method of Keyfitz and Flieger (1990).
#'
#' @param b0 Two-element vector of regression coefficients provided in Keyfitz and Flieger (1990).
#' @param nMx Period central death rates = nDx/nKx.
#'
#' @note
#' Utility function called by \code{\link{life.Leslie}}
#'
#' @return
#' The first two values (age classes 0-1 and 1-5) of the nax column of a period life table.
#'
#' @references
#' Keyfitz, N. and Flieger, W. (1990). World population growth and aging: Demographic trends in the late twentieth century. Chicago: University of Chicago Press.
#'
#' @seealso \code{\link{life.Leslie}}.
#'
#' @keywords keyfitz life-table leslie-matrix
#'
#' @export

keyfitz <- function(b0, nMx){
        nax12 <- c(0,0)
        nax12[1] <- b0[1] + b0[2] * nMx[1]
        nax12[2] <- 1.5
        return(nax12)
}


#' @title
#' Life table for Leslie matrix projections.
#'
#' @description
#' Constructs either a period or cohort life table from enumerated deaths and mid-interval population estimates.
#'
#' @param x Age at the beginning of the age classes of the life table.
#' @param nDx Number of deaths.
#' @param nKx Population size.
#' @param b0 Coefficients used in Keyfitz-Flieger graduation.
#' @param b1 First set of coefficients used in Coale-Demeny graduation.
#' @param b4 Second set of coefficients used in Coale-Demeny graduation.
#' @param type Type of life table calculation: "kf", "cd", or "cohort". Default is "kf".
#' @param nxx individuals-years lived by those dying in the last (possibly open) age-class. If nxx=0,
#' the person-years lived by those dying in the interval is the inverse of the central
#' death rate (corresponding to exponentially distributed failure times).
#' @param iwidth Width of the age intervals.
#' @param width12 Width of the first two age classes.
#'
#' @details
#' Constructs a period or cohort life tables from enumerated deaths and mid-interval population sizes
#' (period) or enumerated deaths and person-years at risk (cohort). x, nDx, and nKx must all the be
#' same length.
#'
#' There are currently three options for life table construction. The first two are for the construction
#' of period life tables. They differ only in the way that person-years lived by those dying in the first
#' two intervals are handled. For type="kf", the default, the first two values of nax estimated using
#' Keyfitz and Fleiger’s (1990) regression method. For type="cd", Coale and Demeny’s method
#' (1983) is used. The Coale-Demeny method uses different coefficients depending on the level of
#' early mortality. As a result, this method may work better for high-mortality populations.
#'
#' The third type of life table is a cohort life table, for which the conversion from mortality rates to
#' probabilities of death is unnecessary, so the nax column of the life table is of limited interest.
#'
#' @return
#' A dataframe with nine columns:
#' \item{x}{Age at the beginning of the interval.}
#' \item{nax}{Individuals-years lived by those dying in the interval x to x + n.}
#' \item{nMx}{Period central death rate.}
#' \item{nqx}{Probability of death between ages x and x + n.}
#' \item{lx}{Probability of survival to exact age x.}
#' \item{ndx}{Proportion of deaths occurring between ages x and x + n.}
#' \item{nLx}{Individuals-years lived in the interval x to x + n.}
#' \item{Tx}{Individuals-years of life left in the cohort at age x.}
#' \item{ex}{Life expectancy at age x.}
#'
#' @note
#' Calls functions \code{\link{keyfitz}}, \code{\link{coale}} or \code{\link{cohort}}
#'
#' @references
#' Keyfitz, N. (1977). Introduction to the mathematics of populations. 2nd ed. Menlo Park: Addison-Wesley.
#'
#' Coale, A., Demeny, P. and Vaughn, B. (1983). Regional model life tables and stable populations. 2nd ed. New York: Academic Press.
#'
#' Keyfitz, N. and Flieger, W. (1990). World population growth and aging: Demographic trends in the late twentieth century. Chicago: University of Chicago Press.
#'
#' Preston, S.H. Heuveline, P. and Guillot, F. (2001). Demography: Measuring and modeling population processes. Oxford: Blackwell.
#'
#' @seealso \code{\link{keyfitz}} \code{\link{coale}} \code{\link{cohort}} \code{\link{Leslie.matrix}}
#'
#' @keywords life-table leslie-matrix
#'
#' @export


life.Leslie <- function(x, nDx, nKx, b0 = c(0.07, 1.7), b1 = c(0.053, 2.8),
                        b4 = c(1.522, 1.518), type = "kf", nxx = 0, iwidth = 5, width12 = c(1, 4)){

        nmax <- length(x)

        # fixxing low values
        nKx[nKx <= 0] <- 0.0001
        nDx[nDx <= 0] <- 0.0001
        nDx[nmax] <- nKx[nmax]


        # period rates n_M_x
        nMx <- nDx/nKx

        # n_a_x
        nax <- NULL
        nax12 <- switch(type,
                        cd = coale(b1,b4,nMx),
                        kf = keyfitz(b0,nMx),
                        cohort = cohort(width12)
        )
        nax[1] <- nax12[1]
        nax[2] <- nax12[2]

        # nxx closes out the life table.  For low-mortality populations (e.g.,
        # where there is under-reporting of old-age deaths), it is essential
        # that nxx be specified.
        nax[3:(nmax-1)] <- iwidth/2
        if(nxx == 0)  {nax[nmax] <- 1/nMx[nmax]
        } else {nax[nmax] <- nxx}

        n <- c(width12, rep(iwidth, nmax - 3),999) #width of the intervals
        nqx <- (n*nMx) / (1 + (n-nax)*nMx)
        nqx[nmax] <- 1.0

        # preceding is unnecessary for type="cohort"
        if(type=="cohort")   nqx <- nMx

        # survivorship lx
        lx <- cumprod(c(1,1-nqx))

        # l_{x+n}
        lxpn <- lx[-1]
        ndx <- -diff(lx)

        # person-years lived by survivors
        nLx <- n*lxpn + ndx*nax
        Tx <- rev(cumsum(rev(nLx)))
        ex <- Tx/lx[1:nmax]

        # format for printing
        lt <- data.frame(x, nax = round(nax,4),
                         nMx = round(nMx,4),
                         nqx = round(nqx[1:nmax],4),
                         lx = round(lx[1:nmax],4),
                         ndx = round(ndx,4),
                         nLx = round(nLx,4),
                         Tx = round(Tx,2),
                         ex = round(ex,2) )

        return(lt)
}

