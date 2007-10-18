\name{textxy}
\alias{textxy}
\title{Label points in a plot}
\description{
  Function textxy calls function text in order to add text
  to points in a graph. textxy chooses a different position 
  for the text depending on the quadrant. This tends to 
  produces better readable plots.
}
\usage{
textxy(X, Y, labs, cx = 0.5, dcol = "black", m = c(0, 0))
}
\arguments{
  \item{X}{x coordinates of a set of points}
  \item{Y}{y coordinates of a set of points}
  \item{labs}{labels to be placed next to the points}
  \item{cx}{character expansion factor}
  \item{dcol}{colour for the labels}
  \item{m}{coordinates of the origin of the plot (default (0,0))}
}
\value{
  NULL
}
\references{ Graffelman, J. (2006) A guide to biplot calibration. }
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{text}} }
\examples{
x <- runif(50)
y <- runif(50)
plot(x,y)
textxy(x,y,1:50,m=c(mean(x),mean(y)))
}
\keyword{aplot}

