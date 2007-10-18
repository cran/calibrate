"textxy" <-
function (X, Y, labs, cx = 0.5, dcol = "black", m = c(0, 0)) 
{
    posXposY <- ((X >= m[1]) & ((Y >= m[2])))
    posXnegY <- ((X >= m[1]) & ((Y < m[2])))
    negXposY <- ((X < m[1]) & ((Y >= m[2])))
    negXnegY <- ((X < m[1]) & ((Y < m[2])))
    if (sum(posXposY) > 0) 
        text(X[posXposY], Y[posXposY], labs[posXposY], adj = c(-0.3, 
            -0.3), cex = cx, col = dcol)
    if (sum(posXnegY) > 0) 
        text(X[posXnegY], Y[posXnegY], labs[posXnegY], adj = c(-0.3, 
            1.3), cex = cx, col = dcol)
    if (sum(negXposY) > 0) 
        text(X[negXposY], Y[negXposY], labs[negXposY], adj = c(1.3, 
            -0.3), cex = cx, col = dcol)
    if (sum(negXnegY) > 0) 
        text(X[negXnegY], Y[negXnegY], labs[negXnegY], adj = c(1.3, 
            1.3), cex = cx, col = dcol)
}
