bplot <- function (Fr, G, rowlab = rownames(Fr), collab = rownames(G), 
          qlt = rep(1, nrow(Fr)), refaxis = TRUE, ahead = T, xl = NULL, 
          yl = NULL, frame = F, qltlim = 0, rowch = 19, colch = 19, 
          qltvar = NULL, rowcolor = "red", colcolor = "blue", 
          rowmark = TRUE, colmark = TRUE, rowarrow = FALSE, colarrow = TRUE, 
          markrowlab = TRUE, markcollab = TRUE, xlab = "", ylab = "", 
          cex.rowlab = 1, cex.rowdot = 0.75, cex.collab = 1, cex.coldot = 0.75, 
          cex.axis = 0.75, lwd = 1, arrowangle = 10, ...) 
{
  opar <- par()
  if (is.null(xl)) 
    xl <- c(min(Fr[, 1]) - (max(Fr[, 1]) - min(Fr[, 1]))/20, 
            max(Fr[, 1]) + (max(Fr[, 1]) - min(Fr[, 1]))/20)
  if (is.null(yl)) 
    yl <- c(min(Fr[, 2]) - (max(Fr[, 2]) - min(Fr[, 2]))/20, 
            max(Fr[, 2]) + (max(Fr[, 2]) - min(Fr[, 2]))/20)
  if (frame == F) {
    opar <- par(yaxt = "n", xaxt = "n", bty = "n", 
                pch = 16)
  }
  plot(0, 0, asp = 1, type = "n", xlim = xl, ylim = yl, 
       cex.lab = cex.axis, cex.axis = cex.axis, xlab = "", 
       ylab = "", ...)
  title(xlab = xlab, line = 2, cex.lab = cex.axis)
  title(ylab = ylab, line = 2, cex.lab = cex.axis)
  if (refaxis) {
    abline(h = 0)
    abline(v = 0)
  }
  if (rowmark) {
    points(Fr[, 1], Fr[, 2], pch = rowch, col = rowcolor, 
           cex = cex.rowdot * qlt, asp = 1)
  }
  if (colmark) {
    if (is.null(qltvar)) 
      points(G[, 1], G[, 2], cex = cex.coldot, pch = colch, 
             col = colcolor)
    else points(G[qltvar > qltlim, 1], G[qltvar > qltlim, 
                                         2], cex = cex.coldot, pch = colch, col = colcolor, 
                asp = 1)
  }
  if (colarrow) {
    if (ahead) {
      arrows(0, 0, G[, 1], G[, 2], length = 0.1, col = colcolor, 
             lwd = lwd, angle=arrowangle)
    }
    else {
      if (is.null(qltvar)) 
        arrows(0, 0, G[, 1], G[, 2], length = 0, lwd = lwd, angle=arrowangle)
      else arrows(0, 0, G[qltvar > qltlim, 1], G[qltvar > 
                                                   qltlim, 2], 
                  length = 0.1, col = colcolor, lwd = lwd, angle = arrowangle)
    }
  }
  if (rowarrow) {
    if (ahead) {
      arrows(0, 0, Fr[, 1], Fr[, 2], length = 0.1, col = rowcolor, 
             lwd = lwd, angle = arrowangle)
    }
    else {
      if (is.null(qltvar)) 
        arrows(0, 0, Fr[, 1], Fr[, 2], length = 0, lwd = lwd, angle = arrowangle)
      else arrows(0, 0, Fr[qltvar > qltlim, 1], Fr[qltvar > 
                                                     qltlim, 2], 
                  length = 0.1, col = rowcolor, lwd = lwd, angle = arrowangle)
    }
  }
  if (markrowlab) {
    textxy(Fr[, 1], Fr[, 2], rowlab, cex = cex.rowlab)
  }
  if (markcollab) {
    if (is.null(qltvar)) 
      textxy(G[, 1], G[, 2], collab, cex = cex.collab)
    else textxy(G[qltvar > qltlim, 1], G[qltvar > qltlim, 
                                         2], collab[qltvar > qltlim], cex = cex.collab)
  }
  par(opar)
}
