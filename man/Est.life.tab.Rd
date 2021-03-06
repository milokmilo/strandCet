% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Est.life.tab.R
\name{Est.life.tab}
\alias{Est.life.tab}
\title{Estimated life table.}
\usage{
Est.life.tab(Est.qx, age = age, n = 1000)
}
\arguments{
\item{Est.qx}{Estimated mortality vector with same length as age.}

\item{age}{The age at the beginning of the age classes of the life table.}

\item{n}{The amount of individuals of the first age class in the theoretical population.}
}
\value{
A dataframe with seven columns:
\item{age}{Age at the beginning of the interval.}
\item{qx}{Probability of death between ages x and x + n.}
\item{nx}{Number of survivors at age x in a theoretical cohort starting with n individuals.}
\item{dx}{Number of deaths at age x in a theoretical cohort starting with n individuals.}
\item{lx}{Probability of survival to exact age x.}
\item{ex}{Life expectancy at age x.}
\item{Zx}{Instant death rate at age x}
}
\description{
Performs a regular life table from estimated and smoothed mortality rates at age.
}
\details{
Constructs a cohort life table from estimated and smoothed mortality rates calculated with mortality
functions as \code{\link{Si.mod}} or \code{\link{HP.mod}}
}
\examples{

modSi <- Si.mod(data = cetaceans, rm = 2,
                par = c(0.3159462,  0.1860541, -1.2802880,  1.1733226,  0.0170314))

dataSi <- Si.pred(data = cetaceans, Sout = modSi, rm = 2)

Est.life.tab(Est.qx = dataSi$qx.tot, age = 0:29, n = 1000)

}
\references{
Preston, S.H., Heuveline, P. and Guillot, F. (2001). Demography: Measuring and modeling population processes. Oxford: Blackwell.
}
\seealso{
\code{\link{Si.mod}} \code{\link{HP.mod}}
}
\keyword{Heligman-Pollard}
\keyword{Siler}
\keyword{life-table}
