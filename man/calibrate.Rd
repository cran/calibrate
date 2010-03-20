\name{calibrate}
\alias{calibrate}
\title{Calibration of Biplot and Scatterplot Axis}
\description{Routine for the calibration of any axis (variable vector) in a biplot or a scatterplot}
\usage{
calibrate(g,y,tm,Fr,tmlab=tm,tl=0.05,dt=TRUE,dp=FALSE,
    lm=TRUE,verb=TRUE,axislab="",reverse=FALSE,
    alpha=NULL,labpos=1,weights=diag(rep(1,length(y))),
    axiscol="blue",cex.axislab=0.75,graphics=TRUE,where=3,
    laboffset=c(0,0),m=matrix(c(0,0),nrow=1),markerpos=3,
    showlabel=TRUE,lwd=1,shiftvec=c(0,0),shiftdir="none",shiftfactor=1.05) 
}
\arguments{
  \item{g}{ the vector to be calibrated (2 x 1). }
  \item{y}{ the data vector corresponding to \code{g}, appropriately centred and/or standardized. }
  \item{tm}{ the vector of tick marks, appropiately centred and/or scaled. }
  \item{Fr}{ the coordinates of the rows markers in the biplot. }
  \item{tmlab}{ a list or vector of tick mark labels. }
  \item{tl}{ the tick length. By default, the tick markers have length 0.05.}
  \item{dt}{ draw ticks. By default, ticks markers are drawn. Set dt=F in order to compute
calibration results without actually drawing the calibrated scale. }
  \item{dp}{ drop perpendiculars. With dp=T perpendicular lines will be drawn from the row markers 
  specified by Fr onto the calibrated axis. This is a graphical aid to read off the values in the
  corresponding scale. }
  \item{lm}{ label markers. By default, all tick marks are labelled. Setting lm=F turns off the labelling
of the tick marks. This allows for creating tick marks without labels. It is particularly useful for
creating finer scales of tickmarks without labels.}
  \item{verb}{ verbose parameter (F=be quiet, T=show results). }
  \item{axislab}{ a label for the calibrated axis. }
  \item{reverse}{ puts the tick marks and tick mark labels on the other side of the axis. }
  \item{alpha}{ a value for the calibration factor. This parameter should only be specified if
 a calibration is required that is different from the one that is optimal for data recovery.}
  \item{labpos}{ position of the label for the calibrated axis (1,2,3 or 4). }
  \item{laboffset}{ offset vector for the axis label. If specified, shifts the label by the specified amounts with respect to the current position. }
  \item{weights}{ a matrix of weights (optional).  }
  \item{axiscol}{ color of the calibrated axis. }
  \item{cex.axislab}{ character expansion factor for axis label and tick mark labels. }
  \item{graphics}{ do graphics or not (F=no graphical output, T=draws calibrated scale). }
  \item{where}{ label placement (1=beginning,2=middle,3=end). }
  \item{m}{ vector of means. }
  \item{markerpos}{ position specifier for the tick mark labels (1,2,3 or 4). }
  \item{showlabel}{ show axis label in graph (T) or not (F). }
  \item{lwd}{ line with for the calibrated axis }
  \item{shiftvec}{ a shift vector for the calibrated axis ((0,0) by
    default) }
  \item{shiftdir}{ indicates in which direction the axis should be
    shifted ("left","right" or "none"). This direction is w.r.t. vector \code{g}}
  \item{shiftfactor}{ scalar by which the shift vector is stretched (or
    shrunken). By default, the length of the shift vector is stretched
    by 5 percent (\code{shiftfactor} = 1.05)}
}
\details{
  This program calibrates variable vectors in biplots and scatterplots, by drawing tick marks along 
  a given the vector and labelling the tick marks with specified values. The optimal calibration is
  found by (generalized) least squares. Non-optimal calibrations are possible by specifying a 
  calibration factor (alpha).
}
\value{
  Returns a list with calibration results
  \item{useralpha}{calibration factor specified by the user}
  \item{optalpha}{optimal calibration factor}
  \item{lengthoneunit}{length in the plot of one unit in the scale of the calibrated variable}
  \item{gof}{goodness of fit (as in regression)}
  \item{gos}{goodness of scale}
  \item{M}{coordinates of the tick markers}
  \item{ang}{angle in degrees of the biplot axis with the positive
    x-axis}
  \item{shiftvec}{the supplied or computed shift vector}
  \item{yt}{fitted values for the variable according to the calibration}
  \item{e}{errors according to the calibration}
  \item{Fpr}{coordinates of the projections of the row markers onto the calibrated axis}
  \item{Mn}{coordinates of the tick marker end points}
}
\references{
Gower, J.C. and Hand, D.J., (1996) Biplots. Chapman & Hall, London

Graffelman, J. and van Eeuwijk, F.A. (2005) Calibration of multivariate scatter plots for 
exploratory analysis of relations within and between sets of variables in genomic research
Biometrical Journal, 47(6) pp. 863-879.

Graffelman, J. (2006) A guide to biplot calibration.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{biplot}}}
\examples{
x <- rnorm(20,1)
y <- rnorm(20,1)
x <- x - mean(x)
y <- y - mean(y)
z <- x + y
b <- c(1,1)
plot(x,y,asp=1,pch=19)
tm<-seq(-2,2,by=0.5)
Calibrate.z <- calibrate(b,z,tm,cbind(x,y),axislab="Z",graphics=TRUE)
}
\keyword{multivariate}
