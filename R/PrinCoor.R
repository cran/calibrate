PrinCoor <- function(Dis, eps = 1e-10) 
{
  n <- nrow(Dis)
  A <- -0.5 * (Dis * Dis)
  n <- nrow(Dis)
  I <- diag(rep(1, n))
  H <- I - (1/n) * ones(n, 1) %*% ones(1, n)
  B <- H %*% A %*% H
  Out <- eigen(B)
  Dl <- diag(Out$values)
  la <- diag(Dl)
  la2 <- la*la
  k <- sum(la > eps)
  
  lar <- round(la,digits=10)
  npos <- sum(lar>0)
  nneg <- sum(lar<0)
  nzer <- sum(lar==0)
  cat("There are",npos,"positive, ",nneg,"negative, ",nzer,"zero eigenvalues\n")
  
  noverlargestnegative <- sum(la >= abs(la[n]))
  cat("There are",noverlargestnegative,
      "positive eigenvalues exceeding the largest negative\n")
  
  Dk <- Dl[1:k, 1:k]
  V <- Out$vectors
  Vk <- V[, 1:k]
  if (k == 1) {
    X <- Vk * sqrt(Dk)
    X <- matrix(X, ncol = 1)
  }
  else X <- Vk %*% sqrt(Dk)
  
  rownames(X) <- rownames(Dis)
  Y <- V%*%Dl
  rownames(Y) <- rownames(Dis)
  
  Qd <- X*X
  Qb <- Y*Y
  
  #
  # Overall goodness-of-fit
  #
  
  #
  # standard
  #
  
  total <- sum(la)
  fr <- la/total
  cu <- cumsum(fr)
  standard.decom <- cbind(la, fr, cu)
  colnames(standard.decom) <- c("Lambda","Fraction","Cumulative")
  
  #
  # positive only
  #
  
  indpos <- la > 0
  posev <- la[indpos]
  posfr <- posev/sum(posev)
  poscu <- cumsum(posfr)
  positive.decom <- cbind(posev, posfr, poscu)
  colnames(positive.decom) <- c("Lambda","Fraction","Cumulative")
  
  #
  # absolute values
  #
  
  absev <- abs(la)
  absfr <- absev/sum(absev)
  abscu <- cumsum(absfr)
  absolute.decom <- cbind(absev, absfr, abscu)
  colnames(absolute.decom) <- c("Lambda","Fraction","Cumulative")
  
  #
  # squared eigenvalues
  #
  
  squaredfr <- la2/sum(la2)
  squaredcu <- cumsum(squaredfr)
  squared.decom <- cbind(la2, squaredfr, squaredcu)
  colnames(squared.decom) <- c("Lambda-squared","Fraction","Cumulative")
  
  #
  # point and pairwise goodness of fit calculations; Euclidean case
  #
  
  w_i <- rowSums(Qd)
  g_i <- (Qd[,1]+Qd[,2])/w_i # equation 9 
    
  D2 <- as.matrix(dist(X[, 1:2]))
  diag(D2) <- 1
  
  Dcopy <- Dis
  diag(Dcopy) <- 1
  
  g_ij <- D2/Dcopy # equation 11
  
  #
  # Euclidean subspace
  #
  
  wl_i <- rowSums(Qd[,1:noverlargestnegative])
  gl_i <- rowSums(Qd[,1:2])/wl_i # equation 12
  
  
  Dell <- as.matrix(dist(X[, 1:noverlargestnegative]))
  diag(Dell) <- 1
  gl_ij <- D2/Dell # equation 13
  
  w.squared_i <- rowSums(Qb)
  gb_i <- rowSums(Qb[,1:2])/w.squared_i # equation 17

  #
  # Error statistics
  #
  
  Er <- Dis-D2
  Er2 <- Er*Er
  total.ess <- sum(Er2[lower.tri(Er2)])
  
  gt_ij <- Er2/total.ess # equation 18
  
  gt_i <- apply(Er2,1,sum)
  toti <- sum(gt_i)
  gt_i <- gt_i/toti  
  
  rn <- rownames(Dis)
  npairs <- 0.5*n*(n-1)
  PairStats <- matrix(NA,nrow=npairs,ncol=6)
  PairStats <- data.frame(PairStats)
  counter <- 1
  for(i in 2:n) {
    for(j in 1:(i-1)) {
      PairStats[counter,1] <- rn[i]
      PairStats[counter,2] <- rn[j]
      PairStats[counter,3] <- g_ij[i,j] # 11
      PairStats[counter,4] <- gl_ij[i,j] # 13
      PairStats[counter,5] <- Er2[i,j] # 18
      counter <- counter + 1
    }
  }
  
  tess <- sum(PairStats[,5]) # total error sum-of-squares
  e.contr.pair <- PairStats[,5]/tess
  
  PairStats[,6] <- e.contr.pair
  
  colnames(PairStats) <- c("ID1","ID2","g_ij","gl_ij",
                              "err2","contr_err")
  
  
  RowStats <- data.frame(g_i,gl_i,gb_i,gt_i)
  
  return(list(X = X, la = la, B = B, standard.decom = standard.decom, 
              positive.decom = positive.decom, 
              absolute.decom = absolute.decom, 
              squared.decom = squared.decom,
              RowStats = RowStats,
              PairStats = PairStats))
}
