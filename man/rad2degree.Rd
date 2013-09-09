\name{rad2degree}
\alias{rad2degree}
\title{
   Convert radians to degrees.
}
\description{
\code{rad2degree converts radians to degrees.}
}
\usage{
rad2degree(x)
}
\arguments{
  \item{x}{an angle in radians}
}
\value{
  the angle with the positive x-axis in degrees.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
   x <- pi/2
   a <- rad2degree(x)
   cat("angle is",a,"degrees\n")
}
\keyword{arith}

