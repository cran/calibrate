\name{ones}
\alias{ones}
\title{Generates a matrix of ones}
\description{
   \code{ones} generates a matrix of ones. 
}
\usage{
ones(n, p = n)
}
\arguments{
  \item{n}{ number of rows }
  \item{p}{ number of columns }
}
\details{
  if only n is specified, the resulting matrix will be square.
}
\value{
  a matrix filled with ones.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{matrix}}}
\examples{
Id <- ones(3)
print(Id)
}
\keyword{multivariate}

