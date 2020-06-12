\name{bplot}
\alias{bplot}
\title{General function for making biplots }
\description{
  Function \code{bplot} creates biplots on the basis matrices of row
  and column markers.
}
\usage{
bplot(Fr,G,rowlab=rownames(Fr),collab=rownames(G),qlt=rep(1,nrow(Fr)),
refaxis=TRUE,ahead=T,xl=NULL,yl=NULL,frame=F,qltlim=0,rowch=19,
colch=19,qltvar=NULL,rowcolor="red",colcolor="blue",rowmark=TRUE,
colmark=TRUE,rowarrow=FALSE,colarrow=TRUE,markrowlab=TRUE,
markcollab=TRUE,xlab="",ylab="",cex.rowlab=1,cex.rowdot=0.75,
cex.collab=1,cex.coldot=0.75,cex.axis=0.75,lwd=1,arrowangle=10,...)}
\arguments{
  \item{Fr}{matrix with coordinates of the row markers.}
  \item{G}{matrix with coordinates of the column markers.}
  \item{rowlab}{vector with labels for the rows.}
  \item{collab}{vector with labels for the columns.}
  \item{qlt}{goodness of fit of the rows.}
  \item{refaxis}{draw coordinate system \code{refaxis=TRUE} or not.}
  \item{ahead}{put a head on the vectors \code{ahead=TRUE} or not.}
  \item{xl}{limits for the x-axis.}
  \item{yl}{limits for the y-axis.}
  \item{frame}{draw a box around the plot \code{frame=TRUE} or not.}
  \item{qltlim}{draw only the vectors with a goodness of fit larger than \code{qltlim}.}
  \item{rowch}{character used for the row markers.}
  \item{colch}{character used for the column markers.}
  \item{qltvar}{vector with the goodness of fit of each variable.}
  \item{rowcolor}{colour used for the row markers.}
  \item{colcolor}{colour used for the column markers.}
  \item{rowmark}{show row markers (\code{rowmark=TRUE}) or not.}
  \item{colmark}{show column markers (\code{colmark=TRUE}) or not.}
  \item{rowarrow}{draw vectors from the origin to the row markers (\code{rowarrow=TRUE}) or not.}
  \item{colarrow}{draw vectors from the origin to the column markers (\code{colarrow=TRUE}) or not.}
  \item{markrowlab}{depict row marker labels (\code{rowlab=TRUE}) or not.}
  \item{markcollab}{depict column marker labels (\code{collab=TRUE}) or not.}
  \item{xlab}{a label for the x-axis.}
  \item{ylab}{a label for the y-axis.}
  \item{cex.rowlab}{expansion factor for the row labels.}
  \item{cex.rowdot}{expansion factor for the row markers.}
  \item{cex.collab}{expansion factor for the column labels.}
  \item{cex.coldot}{expansion factor for the column markers.}
  \item{cex.axis}{expansion factor for the axis.}
  \item{lwd}{line width for biplot vectors.}
  \item{arrowangle}{angle for the edges of the arrowhead.}
  \item{\dots}{extra arguments for plot.}
}
\value{None. The function produces a graphic.}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
set.seed(123)
X <- matrix(runif(40),byrow=TRUE,ncol=4)
colnames(X) <- paste("X",1:ncol(X),sep="")
out.pca <- princomp(X,cor=TRUE)
Fp <- out.pca$scores
Gs <- as.matrix(unclass(out.pca$loadings))
bplot(Fp,Gs,colch=NA)
}
\keyword{multivariate}

