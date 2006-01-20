\name{dlines}
\alias{dlines}
\title{Connect two sets of points by lines}
\description{
  \code{dlines} connects two sets of points by lines in
  a rowwise manner.
}
\usage{
dlines(SetA, SetB, lin = "dotted")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SetA}{matrix with the first set of points}
  \item{SetB}{matrix with teh second set of points}
  \item{lin}{linestyle for the connecting lines}
}
\value{
  NULL
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{lines}}}
\examples{
X <- matrix(runif(20),ncol=2)
Y <- matrix(runif(20),ncol=2)
plot(rbind(X,Y))
text(X[,1],X[,2],paste("X",1:10,sep=""))
text(Y[,1],Y[,2],paste("Y",1:10,sep=""))
dlines(X,Y)
}
\keyword{aplot}

