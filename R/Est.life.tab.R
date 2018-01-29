#' @title
#' Estimated life table.
#'
#' @description
#' Performs a regular life table from estimated and smoothed mortality rates at age.
#'
#' @param Est.qx Estimated mortality vector with same length as age.
#' @param age The age at the beginning of the age classes of the life table.
#' @param n The amount of individuals of the first age class in the theoretical population.
#'
#' @details
#' Constructs a cohort life table from estimated and smoothed mortality rates calculated with mortality
#' functions as \code{\link{Si.mod}} or \code{\link{HP.mod}}
#'
#' @return
#' A dataframe with seven columns:
#' \item{age}{Age at the beginning of the interval.}
#' \item{qx}{Probability of death between ages x and x + n.}
#' \item{nx}{Number of survivors at age x in a theoretical cohort starting with n individuals.}
#' \item{dx}{Number of deaths at age x in a theoretical cohort starting with n individuals.}
#' \item{lx}{Probability of survival to exact age x.}
#' \item{ex}{Life expectancy at age x.}
#' \item{Zx}{Instant death rate at age x}
#'
#' @references
#' Preston, S.H., Heuveline, P. and Guillot, F. (2001). Demography: Measuring and modeling population processes. Oxford: Blackwell.
#'
#' @seealso \code{\link{Si.mod}} \code{\link{HP.mod}}
#'
#' @keywords life-table Siler Heligman-Pollard
#'
#' @examples
#'
#' modSi <- Si.mod(data = cetaceans, rm = 2,
#'                 par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))
#'
#' dataSi <- Si.pred(data = cetaceans, Sout = modSi, rm = 2)
#'
#' Est.life.tab(Est.qx = dataSi$qx.tot, age = 0:29, n = 1000)
#'
#' @export

Est.life.tab <- function (Est.qx, age = age, n = 1000){

        lifetab <- data.frame(age = age, qx = Est.qx)
        n <- n

        # nx and dx
        for (j in 1:nrow(lifetab)){
                if (j == 1) {dx <- lifetab$qx[1] * n; nx <- n
                } else {
                        nx <- c(nx, nx[j - 1] - dx[j - 1])
                        dx <- c(dx, lifetab$qx[j] * nx[j])
                }
        }
        lifetab <- data.frame(lifetab, nx = nx, dx = dx)

        # lx - Survivorship-at-age percent
        lifetab$lx <- lifetab$nx/n

        # ex - Life expentancy at age ex = Sumlx/lx
        ex <- vector("numeric")
        for (j in 1:nrow(lifetab)) {
                e <- sum(lifetab$lx[j:nrow(lifetab)])/lifetab$lx[j]
                ex <- c(ex, e) }
        lifetab$ex <- ex

        # Zx - Total mortality-at-age -L(nt/no)/t
        Zx <- c(NA)
        for (j in 1:nrow(lifetab)){
                if (j == 1) {Zx <- vector("numeric")}
                z <- -log(lifetab$nx[j + 1]/lifetab$nx[j])/1
                # Correction for the last mortality
                if (j == nrow(lifetab)) {z <- 1}
                Zx <- c(Zx,z)
        }
        lifetab$Zx <- Zx

        # Format for printing
        lifetab$qx <- round(lifetab$qx,3)
        lifetab$nx <- round(lifetab$nx,3)
        lifetab$dx <- round(lifetab$dx,3)
        lifetab$lx <- round(lifetab$lx,3)
        lifetab$ex <- round(lifetab$ex,3)
        lifetab$Zx <- round(lifetab$Zx,3)

        return(lifetab)
}

