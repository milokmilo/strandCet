#' @title
#' Leslie matrix
#'
#' @description
#' Generates a Leslie matrix for demographic projection from vectors of age-specific cumulative survival and fertility.
#'
#' @param lx A vector of either age-specific cumulative survival or person-years lived in the interval
#' @param mx Age-specific fertility rates.
#' @param L Logical. If ’TRUE’, lx is taken to be individuals-years lived in the interval nLx,
#' while if ’FALSE’, lx is taken to be cumulative survival to exact age x + n.
#' @param peryear Multiplier for fertility.
#' @param one.sex Logical. If ’TRUE’, fertility rates will be divided by 1/(1+SRB).
#' @param SRB Sex ratio at birth.
#' @param infant.class Logical. ’TRUE’ if lx contains a value for the infant age-class.
#'
#' @details
#' Constructs a k x k age-classified demographic projection matrix with age-specific survival probabilities
#' along the sub-diagonal and age-specific fertilities along the first row of the matrix.
#'
#' lx and mx are assumed to be of the same length. The resulting matrix is truncated to insure that
#' there are no post-reproductive classes. This is important for ensuring irreducibility of the resulting
#' matrix.
#'
#' If mx is longer than lx, mx is trucated to be the same length as lx. If lx is longer than mx, a warning
#' is issed and lx is truncated to be the same length as mx.
#'
#' Fertility is assumed to be birth-flow (Caswell 2001). That is, breeding is assumed to be continuous
#' and the individual elements of the first row of the Leslie matrix are averaged over successive ageclasses.
#' Fertility rates are typically given in annualized form. If this is the case and the age-classes
#' are wider than one year, then peryear can be used to appropriately scale up the annual values.
#'
#' The default behavior is to use person-years lived in the interval as the survival measure. If infant.class=TRUE,
#' lx is taken to have a value for the infant age class (i.e., a shorter class width than the other elements of lx.
#' What is done when there is an infant class depends on what the values in lx represent.
#' If L=TRUE, then the first two values of lx are combined to form the total person-years for the first ageclass
#' in the Leslie matrix. Human demographic data from abridged life tables typically come with
#' age classes x = 0, 1, 5, 10, ... Thus, combining the person-years for the first two age classes gives
#' an initial age class of the correct width. If infant.class=TRUE and L=FALSE, the second element
#' of lx is deleted. Creating a Leslie matrix from other forms of non-standard early age-classes can
#' be accomplished by pre-processing lx and using the option infant.class=FALSE.
#'
#' The human sex ratio at birth (male births/female births) is remarkably close to SRB=1.05 across a
#' wide range of populations and this is the default value for SRB.
#'
#' The resulting matrix has class "leslie.matrix". This class is not used extensively but will be in future
#' evelopment.
#'
#' @return
#' A k x k age-classified demographic projection matrix with class "leslie.matrix".
#'
#' @references
#' Keyfitz, N. (1977). Introduction to the mathematics of populations. 2nd ed. Menlo Park: Addison-Wesley.
#'
#' Preston, S.H., Heuveline, P. and Guillot, F. (2001). Demography: Measuring and modeling population processes. Oxford: Blackwell.
#'
#' Caswell, H. (2001). Matrix population models: Construction, analysis, and interpretation. 2nd ed. Sunderland, MA: Sinauer.
#'
#' @keywords Leslie-matrix
#'
#' @export

Leslie.matrix <- function (lx, mx, L = TRUE, peryear = 5, one.sex = TRUE, SRB = 1,
                           infant.class = TRUE)
{
        lx[lx <= 0] <- 0.0001

        len1 <- length(lx)
        len2 <- length(mx)
        if (len1 > len2) {
                warning("length of lx greater than the length of mx,\n lx truncated to length of mx")
                lx <- lx[1:len2]
        }
        if (len2 > len1) {
                mx <- mx[1:len1]
        }
        if (infant.class){
                mx <- mx[-2]
        }
        fages <- which(mx > 0)
        k <- max(fages)
        mx <- mx[1:k]

        if (L) {
                L1 <- lx[1] + lx[2]
                s <- L1
                lx <- c(L1, lx[-c(1, 2)])
        }
        s <- lx[1]
        #  else {
        #    lx <- lx[-2]
        #    s <- sqrt(lx[2]) * peryear
        #  }
        px <- exp(diff(log(lx)))
        px <- px[1:(k - 1)]
        Fx <- NULL
        for (i in 1:k - 1) {
                Fx[i] <- s * (mx[i] + px[i] * mx[i + 1])/2
        }
        Fx <- c(Fx, mx[k])
        if (one.sex)
                Fx <- Fx/(1 + SRB)
        A <- matrix(0, nrow = k, ncol = k)
        A[row(A) == col(A) + 1] <- px
        A[1, ] <- Fx
        class(A) <- "leslie.matrix"
        A
}
