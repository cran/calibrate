\name{origin}
\alias{origin}
\title{Origin}
\description{
  Draws coordinate axes in a plot.
}
\usage{origin(m=c(0,0))}
\arguments{
  \item{m}{ the coordinates of the means (2 x 1). }
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
   X <- matrix(runif(40),ncol=2)
   plot(X[,1],X[,2])
   origin(m=c(mean(X[,1]),mean(X[,2])))
}
\keyword{multivariate}
