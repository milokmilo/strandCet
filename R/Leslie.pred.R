#' @title
#' Project Leslie matrix
#'
#' @description
#' Projects an population vector tmax intervals by pre-multiplication with a Leslie matrix.
#'
#' @param A A k * k projection matrix.
#' @param no A k * 1 population vector.
#' @param tmax Number of time steps to project the vector.
#' @param pop.sum Logical: If ’TRUE’, the age-classes of the projected population are summed, yielding a single total population vector
#'
#' @details
#' Takes an initial population vector, no, and pre-multiplies by the demographic projection matrix,
#' A, tmax times. This projection will be tmax*n years into the future, where n is the width of the
#' age-classes in the Leslie matrix, A.
#'
#' @return
#' If pop.sum=FALSE (the default), the value will be a k x tmax+1 matrix. The first column of the
#' matrix is no and each subsequent column represents the population structure at time step 1, 2, ..., tmax.
#' If pop.sum=TRUE, the value will be a vector of length tmax+1, where each element of the vector is
#' the total population at time t=0, 1, ..., tmax.
#'
#' @references
#' Caswell, H. (2001). Matrix population models: Construction, analysis, and interpretation. 2nd ed. Sunderland, MA: Sinauer.
#'
#' van Groenendael, J., De Kroon, H., Kalisz, S. and Tuljapurkar. S. (1994). Loop analysis: Evaluating life history pathways in population projection matrices. Ecology 75 (8):2410-2415.
#'
#' @seealso \code{\link{Leslie.matrix}} \code{\link{eigen.analysis}}
#'
#' @keywords Leslie-matrix projection
#'
#' @export

Leslie.pred <- function(A, no, tmax = 100, pop.sum = FALSE){
        if(length(no) != dim(A)[1])
                stop("Projection matrix and population vector have different number of states!")
        N <- matrix(0, nrow = length(no), ncol = tmax + 1)
        N[,1] <- no
        pop <- no
        for(t in 1:tmax){
                pop <- A %*% pop
                N[, t + 1] <- pop
        }
        if(pop.sum){
                N <- apply(N, 2, sum)
        }
        return(N)
}

