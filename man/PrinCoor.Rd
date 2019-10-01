\name{PrinCoor}
\alias{PrinCoor}
\title{Function for Principal Coordinate Analysis
}
\description{
  Function \code{PrinCoor} implements Principal Coordinate Analysis, also known as classical metric multidimensional scaling or
  classical scaling. In comparison with other software, it offers refined statistics for goodness-of-fit at the level of individual observations and pairs of observartions.
}
\usage{
PrinCoor(Dis, eps = 1e-10)
}
\arguments{
  \item{Dis}{A distance matrix or dissimilarity matrix}
  \item{eps}{A tolerance criterion for deciding if eigenvalues are zero or not}
}
\details{
 Calculations are based on the spectral decomposition of the scalar product matrix B, derived from the distance matrix.
}
\value{
  \item{X}{The coordinates of the the solution}
  \item{la}{The eigenvalues of the solution}
  \item{B}{The scalar product matrix}
  \item{standard.decom}{Standard overall goodness-of-fit table using all eigenvalues}
  \item{positive.decom}{Overall goodness-of-fit table using only positive eigenvalues}
  \item{absolute.decom}{Overall goodness-of-fit table using absolute values of eigenvalues}
  \item{squared.decom}{Overall goodness-of-fit table using squared eigenvalues}
  \item{RowStats}{Detailed goodness-of-fit statistics for each row}
  \item{PairStats}{Detailed goodness-of-fit statistics for each pair}
}
\references{
Graffelman, J. (2019) Goodness-of-fit filtering in classical metric multidimensional scaling with large
datasets. <doi: 10.1101/708339>

Graffelman, J. and van Eeuwijk, F.A. (2005) Calibration of multivariate scatter plots for 
exploratory analysis of relations within and between sets of variables in genomic research
Biometrical Journal, 47(6) pp. 863-879.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }

\seealso{\code{\link{princomp}}}
\examples{
   data(spaindist)
   results <- PrinCoor(as.matrix(spaindist))
}
\keyword{multivariate}
