\name{shiftvector}
\alias{shiftvector}
\title{
   Compute a shift vector for a calibrated axis.
}
\description{
  \code{shiftvector} computes two shift vectors perpendicular to the
  supplied biplot or scatterplot axis \code{g}. The vector norm is
  computed from the two most extreme data points.
}
\usage{
shiftvector(g, X, x = c(1, 0), verbose = FALSE)
}
\arguments{
  \item{g}{a biplot or scatterplot axis}
  \item{X}{a n by 2 matrix of scatterplot or biplot coordinates}
  \item{x}{reference axis, (1,0) by default}
  \item{verbose}{print information or not}
}
\details{
  \code{shiftvector} locates the tow most extreme datapoints in the
  direction perpendicular to axis \code{g}.
}
\value{
  \item{dr}{the right (w.r.t. the direction of \code{g}) shift vector}
  \item{dl}{the left (w.r.t. the direction of \code{g}) shift vector}
}
\references{
Graffelman, J. and van Eeuwijk, F.A. (2005) Calibration of multivariate scatter plots for 
exploratory analysis of relations within and between sets of variables in genomic research
Biometrical Journal, 47(6) pp. 863-879.

Graffelman, J. (2006) A guide to biplot calibration.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{calibrate}}}

\examples{
X <- matrix(rnorm(100),ncol=2)
Xs <- scale(X)

g <- c(1,1)

plot(Xs[,1],Xs[,2],asp=1,pch=19)
textxy(Xs[,1],Xs[,2],1:nrow(X))

arrows(0,0,g[1],g[2])
text(g[1],g[2],"g",pos=1)

out <- shiftvector(g,X,verbose=TRUE)
dr <- out$dr
dl <- out$dl

arrows(0,0,dl[1],dl[2])
text(dl[1],dl[2],"dl",pos=1)

arrows(0,0,dr[1],dr[2])
text(dr[1],dr[2],"dr",pos=1)
}
\keyword{multivariate}
