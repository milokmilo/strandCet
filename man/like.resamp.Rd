% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HP.mod.R
\name{like.resamp}
\alias{like.resamp}
\title{Local Optimums and Covariance from the optimizer step.}
\usage{
like.resamp(K, log.like.0, opt.cov.d, opt.mu.d, d.keep, d = 10,
  theta.dim = 9)
}
\arguments{
\item{K}{Number of iterations at the importance sampling stage.}

\item{log.like.0}{A vector containing the likelihoods for each row of the prior.}

\item{opt.cov.d}{Covariance matrixes for the local optimums.}

\item{opt.mu.d}{A d x 8 matrix containing the local optimums (sets of parameters from the optimizer step).}

\item{d.keep}{Number of local optimums found in the optimizer step.}

\item{d}{A scalar defining the number of optimizer interations.}

\item{theta.dim}{Number of columns in the prior matrix.}
}
\value{
\item{h.mu}{A d.keep * 8 matrix containing the local optimum result.}
\item{h.sig}{An array with (theta.dim * theta.dim * (K + d.keep)) dimensions containing the covariance matrix for each local optimum.}
\item{log.like}{A vector of likelihoods for each row of H.k.}
}
\description{
Defines some necessary arguments for the function \link{final.resamp}. Removes NAs from the opt.mu.d and opt.cov.d matrixes.
}
\note{
Typically for use immediately before running \link{final.resamp} or within the function \link{HP.mod}
}
\references{
Poole, D. and Raftery, A. (2000). Inference for Deterministic Simulation Models: The Bayesian Melding Approach. Journal of the American Statistical Association 95:1244–1255.

Raftery, A. and Bao, L. (2009). Estimating and Projecting Trends in HIV/AIDS Generalized Epidemics Using Incremental Mixture Importance Sampling. Technical Report 560, Department of Statistics, University of Washington.

Sharrow, D.J., Clark, S.J., Collinson, M.A., Kahn, K. and Tollman, S.M. (2013). The Age Pattern of Increases in Mortality Affected by HIV: Bayesian Fit of the Heligman-Pollard Model to Data from the Agincourt HDSS Field Site in Rural Northeast South Africa. Demogr. Res. 29, 1039–1096.
}
\seealso{
\link{final.resamp} \link{HP.mod}
}
\keyword{Heligman-Pollard}
\keyword{bycatch}
\keyword{likelihood}
\keyword{mortality}
\keyword{resampling}
