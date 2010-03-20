\name{rda}
\alias{rda}
\title{Redundancy analysis}
\description{
  \code{rda} performs redundancy analysis and stores extensive output
  in a list object. 
}
\usage{
rda(X, Y, scaling = 1)
}
\arguments{
  \item{X}{a matrix of x variables}
  \item{Y}{a matrix of y variables}
  \item{scaling}{scaling used for x and y variables. 0: x and y only centered. 1: x and y standardized}
}
\details{
  Results are computed by doing a principal component analyis of the
  fitted values of the regression of y on x.

  Plotting the first two columns of Gxs and Gyp, or of Gxp and Gys provides a biplots of 
  the matrix of regression coefficients.

  Plotting the first two columns of Fs and Gp or of Fp and Gs provides a biplot of the
  matrix of fitted values.
}
\value{
Returns a list with the following results
\item{Yh}{ fitted values of the regression of y on x }
\item{B}{ regression coefficients of the regresson of y on x }
\item{decom}{ variance decomposition/goodness of fit of the fitted values AND of
the regression coefficients }
\item{Fs}{ biplot markers of the rows of Yh (standard coordinates) }
\item{Fp}{ biplot markers of the rows of Yh (principal coordinates) }
\item{Gys}{ biplot markers for the y variables (standard coordinates) }
\item{Gyp}{ biplot markers for the y variables (principal coordinates) }
\item{Gxs}{ biplot markers for the x variables (standard coordinates) }
\item{Gxp}{ biplot markers for the x variables (principal coordinates) }
}
\references{
   Van den Wollenberg, A.L. (1977) Redundancy Analysis, an alternative for canonical 
   correlation analysis. Psychometrika 42(2): pp. 207-219.

   Ter Braak, C. J. F. and Looman, C. W. N. (1994) Biplots in Reduced-Rank Regression.
   Biometrical Journal 36(8): pp. 983-1003.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{princomp}},\code{\link{canocor}},\code{\link{biplot}}}
\examples{
X <- matrix(rnorm(75),ncol=3)
Y <- matrix(rnorm(75),ncol=3)
rda.results <- rda(X,Y)
}
\keyword{multivariate}

