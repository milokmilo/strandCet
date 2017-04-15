#' @title
#' Generation time
#'
#' @description
#' Calculates the generation time for an age or stage-classified demographic projection matrix.
#'
#' @param A Demographic projection matrix.
#' @param peryear Width of the age classes.
#'
#' @details
#' Calculates the generation time (T) for an age or stage-classified demographic projection matrix
#' using the identity:
#'
#' Ro = exp(r*T)
#'
#' where Ro is the net reproduction number and r is the intrinsic rate of increase = log(lambda).
#'
#' Generation time is the amount of time that it takes a typical female to produce Ro offspring or,
#' equivalently, the amount of time it takes a population growing with instantaneous rate r to increase
#' by a factor of Ro.
#'
#' @return
#' The generation time implied by the demographic projection matrix.
#'
#' @note
#' Calls function \code{\link{calc.ro}}, which calculates $R_0$ from the fundamental matrix of the Markov transition
#' matrix (Caswell 2001).
#'
#' @references
#' Keyfitz, N., and Caswell. H. (2005). Applied mathematical demography. 3rd ed. New York: Springer.
#'
#' Caswell, H. (2001). Matrix population models: Construction, analysis, and interpretation. 2nd ed. Sunderland, MA: Sinauer.
#'
#' Preston, S.H., Heuveline, P. and Guillot, F. (2001). Demography: Measuring and modeling population processes. Oxford: Blackwell.
#'
#' @seealso \code{\link{calc.ro}} \code{\link{eigen.analysis}}
#'
#' @keywords leslie-matrix generation
#'
#' @export

gen.time <- function(A, peryear = 5){
        ro <- calc.ro(A)
        ea <- eigen.analysis(A)
        T <- peryear*log(ro)/log(ea$lambda)

        return(T)
}

