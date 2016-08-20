TP_NRWRH <- function(Sa, Sc, adjAC) {
  n.A <- nrow(Sa); n.C <- nrow(Sc)
  Sna <- transferWeight(n.A, n.C, adjAC, flags.lambda)
  Snc <- transferWeight(n.C, n.A, t(adjAC), flags.lambda)
  tmp.Sa <- 1 - (1 - Sa) * (1 - Sna)
  tmp.Sc <- 1 - (1 - Sc) * (1 - Snc)
  
  print(paste0("two pass"))
  mt <- nrwrh_avg(tmp.Sa, tmp.Sc, adjAC)    # size: n.C * n.A
  return( t(mt) )                           # size: n.A * C
}

TP_NRWRH.onePass <- function(Sa, Sc, adjAC) {
  n.A <- nrow(Sa); n.C <- nrow(Sc)
  Sna <- transferWeight(n.A, n.C, adjAC, flags.lambda)
  Snc <- transferWeight(n.C, n.A, t(adjAC), flags.lambda)
  tmp.Sa <- 1 - (1 - Sa) * (1 - Sna)
  tmp.Sc <- 1 - (1 - Sc) * (1 - Snc)
  print("one pass.")
  # 单方向游走
  mt <- nrwrh(tmp.Sa, tmp.Sc, adjAC)
  return(mt) # return: C * A
}