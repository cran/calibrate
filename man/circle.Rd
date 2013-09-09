\name{circle}
\alias{circle}
\title{Draw a circle}
\description{
   \code{circle} draws a circle in an existing plot.
}
\usage{
circle(radius,origin)
}
\arguments{
  \item{radius}{ the radius of the circle }
  \item{origin}{ the origin of the circle }
}
\value{
  NULL
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\examples{
   set.seed(123)
   X <- matrix(rnorm(20),ncol=2)
   plot(X[,1],X[,2])
   circle(1,c(0,0))
}
\keyword{aplot}
