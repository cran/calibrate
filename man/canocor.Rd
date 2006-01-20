\name{canocor}
\alias{canocor}
\title{Canonical correlation analysis}
\description{
  \code{canocor} performs canonical correlation analysis on the
  basis of the standardized variables and stores extensive output
  in a list object.
}
\usage{
canocor(X, Y)
}
\arguments{
  \item{X}{ a matrix containing the X variables }
  \item{Y}{ a matrix containing the Y variables }
}
\details{
  \code{canocor} computes the solution by a singular value 
  decomposition of the transformed between set correlation matrix.
}
\value{
Returns a list with the following results
\item{ccor }{ the canonical correlations }
\item{A }{ canonical weights of the x variables } 
\item{B }{ canonical weights of the y variables } 
\item{U }{ canonical x variates } 
\item{V }{ canonical y variates } 
\item{Fs }{ biplot markers for x variables (standard coordinates) }
\item{Gs }{ biplot markers for y variables (standard coordinates) }
\item{Fp }{ biplot markers for x variables (principal coordinates) }
\item{Gp }{ biplot markers for y variables (principal coordinates) }
\item{fitRxy }{ goodness of fit of the between-set correlation matrix }
\item{fitXs }{ adequacy coefficients of x variables } 
\item{fitXp }{ redundancy coefficients of x variables } 
\item{fitYs }{ adequacy coefficients of y variables } 
\item{fitYp }{ redundancy coefficients of y variables} 
}
\references{ 
   Hotelling, H. (1935) The most predictable criterion. Journal of Educational 
   Psychology (26) pp. 139-142.

   Hotelling, H. (1936) Relations between two sets of variates. Biometrika
   (28) pp. 321-377.

   Johnson, R. A. and Wichern, D. W. (2002) Applied Multivariate Statistical Analysis.
   New Jersey: Prentice Hall.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{cancor}}}
\examples{
set.seed(123)
X <- matrix(runif(75),ncol=3)
Y <- matrix(runif(75),ncol=3)
cca.results <- canocor(X,Y)
}
\keyword{multivariate}

