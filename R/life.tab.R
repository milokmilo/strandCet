#' @title
#' Life table
#'
#' @description
#' Performs a regular life table from a dataframe with ages and number of death individuals.
#'
#' @param x A dataframe withtwo columns: Ages and number of deaths.
#' @param n Inicial number of individuals in the theoretical population.
#'
#' @details
#' Constructs a cohort life table from a dataframe with age classes and enumerated deaths.
#'
#' @return
#' A dataframe with nine columns:
#' \item{age}{Age at the beginning of the interval.}
#' \item{Mx}{Number of observed deaths at age x.}
#' \item{Sx}{Number of survivors at age x in the observed cohort of sum(length(M)) individuals.}
#' \item{nx}{Number of survivors at age x in a theoretical cohort starting with n individuals.}
#' \item{dx}{Number of deaths at age x in a theoretical cohort starting with n individuals.}
#' \item{qx}{Probability of death between ages x and x + n.}
#' \item{lx}{Probability of survival to exact age x.}
#' \item{ex}{Life expectancy at age x.}
#' \item{Zx}{Instant death rate at age x}
#'
#' @references
#' Preston, S.H., Heuveline, P. and Guillot, F. (2001). Demography: Measuring and modeling population processes. Oxford: Blackwell.
#'
#' @keywords life-table
#'
#' @examples
#'
#' life.tab(cetaceans)
#'
#' @export

life.tab <- function(x, n = 1000){

        ## x = data frame with ages and number of deaths
        life <- x
        names(life) <- c("age","Mx")
        n <- n

        ## Sx - Creating survivorship column
        d <- sum(life$Mx) # sum of deaths
        for (j in 1:nrow(life)){
                life$Sx[j] <- d-sum(life$Mx[1:j-1])}

        ## nx - Standardizing survivors [(N/SumM) * 1000]
        life$nx <- (life$Sx/d)*n

        ## dx - Dolphins death-at-age [nx - n(x+1)]
        for (j in 1:nrow(life)) {
                if (j == 1) {dx <- vector("numeric")}
                d <- life$nx[j]-life$nx[j+1]
                if (j == nrow(life)) {d <- life$nx[j]-0}
                dx <- c(dx,d)}
        life$dx <- dx

        ## qx - Death-at-age probability [dx / nx]
        for (j in 1:nrow(life)) {
                if (j == 1) {qx <- vector("numeric")}
                q <- life$dx[j]/life$nx[j]
                qx <- c(qx,q) }
        life$qx <- qx

        ## lx - Survivorship-at-age percent [nx / n]
        life$lx <- life$nx/n

        ## ex - Life expentancy at age ex = Sumly/lx
        ex <- vector("numeric")
        for (j in 1:nrow(life)) {
                e <- sum(life$lx[j:nrow(life)])/life$lx[j]
                ex <- c(ex,e) }
        life$ex <- ex

        ## Zx - Total mortality-at-age -L(nt/no)/t
        Zx <- c(NA)
        for (j in 1:nrow(life)){
                if (j == 1) {Zx <- vector("numeric")}
                z <- -log(life$nx[j+1]/life$nx[j])/1
                if (j == nrow(life)) {z <- 1.00}
                Zx <- c(Zx,z)
        }
        life$Zx <- Zx

        ## Format for printing
        life$qx <- round(life$qx,3)
        life$nx <- round(life$nx,3)
        life$dx <- round(life$dx,3)
        life$lx <- round(life$lx,3)
        life$ex <- round(life$ex,3)
        life$Zx <- round(life$Zx,3)

        return(life)
}
