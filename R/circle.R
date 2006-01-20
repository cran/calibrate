"circle" <-
function(radius=1) {
  x<-seq(-radius,radius,by=0.01)
  y<-sqrt(radius^2-x^2);
  lines(x,y)
  lines(x,-y)
  return(NULL)
}

