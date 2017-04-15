#' @title
#' Caclulate net reproduction number from a demographic projection matrix.
#'
#' @description
#' Calculate the net reproduction number ($R_0$) from an age or stage-classified demographic projection matrix.
#'
#' @param A A demographic projection matrix
#' @param N.out Return the fundamental matrix (N) of the Markov chain.
#'
#' @details
#' Calculates the net reproduction number ($R_0$) from an age or stage-classified demographic projection
#' matrix by first decomposing the k x k projection matrix A into two component matrices, T
#' and F. T collects the transitions between life-cycle stages while F collects the fertility transitions.
#' For an age-classified Leslie matrix, T will contain only the sub-diagonal of A and F will contain
#' only the first row of A. The fundamental matrix is given by $N = (I-T)^-1$, where I is a k x k
#' identity matrix. $R_0$ is the leading eigenvalue of the matrix FN.
#'
#' @return
#' If the (default) option N.out=FALSE is used, the net reproduction number is returned as a single
#' value. If N.out = TRUE, the returned value is a list of two items:
#' \item{ro}{Net reproduction number.}
#' \item{N}{Fundamental matrix.}
#'
#' @references
#' Caswell, H. (2001). Matrix population models: construction, analysis, and interpretation, Second edition. Sinauer, Sunderland, Massachusetts, USA.
#'
#' @seealso \code{\link{Leslie.matrix}}
#'
#' @keywords rho leslie-matrix
#'
#' @export

calc.ro <- function(A, N.out = FALSE){
        # Net reproduction number from Leslie matrix
        # assumes age-structured Leslie matrix
        k <- dim(A)[1]
        T <- A
        T[1,] <- 0                     # matrix of transitions
        F <- matrix(0, nrow = k, ncol = k)
        F[1,] <- A[1,]                 # matrix of births
        N <- solve(diag(k) - T)          # fundamental matrix
        ev <- eigen(F %*% N)
        imax <- which(ev$values == max(Re(ev$values)))
        ro <- ev$values[imax]          # same as FN[1,1]

        if(N.out) out <- list(ro, N)
        else out <- ro

        return(out)
}
