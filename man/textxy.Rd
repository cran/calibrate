\name{textxy}
\alias{textxy}
\title{ Nice placement of labels in a plot}
\description{
  Function \code{textxy} calls function \code{text} in order to add text
  to points in a graph. \code{textxy} chooses a different position 
  for the text depending on the quadrant. This tends to 
  produces better readable plots, with labels fanning away from the origin.
}
\usage{
textxy(X, Y, labs, m = c(0, 0), cex = 0.5, offset = 0.8, ...) 
}
\arguments{
  \item{X}{x coordinates of a set of points}
  \item{Y}{y coordinates of a set of points}
  \item{labs}{labels to be placed next to the points}
  \item{m}{coordinates of the origin of the plot (default (0,0))}
  \item{cex}{character expansion factor}
  \item{offset}{controls the distance between the label and the point. A
  value of 0 will plot labels on top of the point. Larger values give
  larger separation between point and label. The default value is 0.8}
  \item{\dots}{additiona arguments for function \code{text}.}
}
\value{
  NULL
}
\references{ Graffelman, J. (2006) A guide to biplot calibration. }
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{text}} }
\examples{
x <- rnorm(50)
y <- rnorm(50)
plot(x,y,asp=1)
textxy(x,y,1:50,m=c(mean(x),mean(y)))
}
\keyword{aplot}
\keyword{misc}


