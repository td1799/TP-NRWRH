source("utils.R", encoding = "UTF-8")

init.feature.vec <- function(adjAB, i.A) {
  A <- rownames(adjAB); B <- colnames(adjAB)
  total.length <- length(A) + length(B)
  total.names <- c(A, B)
  p.0 <- matrix(0, nrow = total.length, ncol=1, dimnames = list(total.names, c(A[i.A])))
  p.0[A[i.A], 1] <- 1.0
  connected.others <- names(which(adjAB[i.A,] == 1))
  if (length(connected.others) == 0) {
    
  } else {
    p.0[A[i.A], 1] <- eta
    p.0[connected.others, 1] <- (1 - eta) / length(connected.others)
  }
  if (sum(p.0) != 1) {
    stop("initial vector Error.")
  }
  return(p.0)
}

init.feature.matrix <- function(adjAB) {
  A <- rownames(adjAB); B <- colnames(adjAB)
  m.0 <- matrix(0, nrow=length(A)+length(B), ncol=length(A), dimnames = list(c(A, B), A))
  n.A <- length(A)
  for (i in 1:n.A) {
    tmp.con <- names(which(adjAB[i, ]==1))
    tmp.n <- length(tmp.con)
    m.0[i, i] <- 1.0
    if (tmp.n == 0) {
      next
    } else {
      m.0[i, i] <- eta  # 自己留eta
      m.0[tmp.con, i] <- (1 - eta) / tmp.n
    }
  }
  if (sum(m.0) != n.A) {
    stop("Initial matrix Error.")
  }
  return(m.0)
}

build.transfer.matrix <- function(Sa, Sb, adjAB) {
  A <- rownames(adjAB); B <- colnames(adjAB); 
  total.length <- length(A) + length(B)
  total.names <- c(A, B)
  M <- matrix(0, nrow=total.length, ncol=total.length, dimnames=list(total.names, total.names))
  
  diag.A <- diag(1.0 / rowSums(Sa)); dimnames(diag.A) <- list(A, A)
  diag.A[is.infinite(diag.A) | is.na(diag.A)] <- 1.0
  M[A, A] <- diag.A %*% Sa
  M[A, A][names(which(rowSums(adjAB)!=0)), ] <- (1 - lambda) * M[A, A][names(which(rowSums(adjAB)!=0)), ]
  diag.AB <- diag(1.0 / rowSums(adjAB)); dimnames(diag.AB) <- list(A, A)
  diag.AB[is.infinite(diag.AB) | is.na(diag.AB)] <- 1.0
  M[A, B] <- lambda * diag.AB %*% adjAB
  
  diag.B <- diag(1.0 / rowSums(Sb)); dimnames(diag.B) <- list(B, B)
  diag.B[is.infinite(diag.B) | is.na(diag.B)] <- 1.0
  M[B, B] <- diag.B %*% Sb
  M[B, B][names(which(colSums(adjAB)!=0)), ] <- (1 - lambda) * M[B, B][names(which(colSums(adjAB)!=0)), ]
  adjBA <- t(adjAB)
  diag.BA <- diag(1.0 / rowSums(adjBA)); dimnames(diag.BA) <- list(B, B)
  diag.BA[is.infinite(diag.BA) | is.na(diag.BA)] <- 1.0
  M[B, A] <- lambda * diag.BA %*% adjBA
  if (!isTRUE(all.equal(unname(rowSums(M)), rep(1, nrow(M)), tolerance=0.0001))) {
    print("Something wrong in M.")
  }
  return(M)
}

start.transfer.vector <- function(M, p.0) {
  turn <- 0
  p.t <- p.0
  while (TRUE) {
    turn <- turn + 1
    p.t1 <- (1 - restart.ratio)* (t(M)%*%p.t) + restart.ratio * p.0
    L1 <- L1_norm(p.t1 - p.t)
    if (L1 < LIMIT_TOR) {
      break
    }
    p.t <- p.t1
  }
  return(p.t1)
}

# random walk progress. for matrix
# parameters:
#   M   : transition matrix, [Mcc, Mcp, Mcd] the first row.
#   m.0  : initial matrix
#   r   : restart ratio
start.transfer.matrix <- function(M, m.0) {
  turn <- 0
  m.t <- m.0
  while (TRUE) {
    turn <- turn + 1
    m.t1 <- (1 - restart.ratio) * (t(M) %*% m.t) + restart.ratio * m.0
    L1 <- L1_norm(m.t1 - m.t)
    if (L1 < LIMIT_TOR) {
      break
    }
    m.t <- m.t1
  }
  return(m.t1)
}

# Sa, Sc should have been already adjusted somewhere else.
nrwrh <- function(Sa, Sc, adjAC) {
  M <- build.transfer.matrix(Sa, Sc, adjAC)
  m0 <- init.feature.matrix(adjAC)
  mt <- start.transfer.matrix(M, m0)
  C <- rownames(Sc)
  return(mt[C,]) # size: n.C * n.A
}

nrwrh_avg <- function(Sa, Sc, adjAC) {
  mt.A2C <- nrwrh(Sa, Sc, adjAC)     # size: n.C * n.A
  mt.C2A <- nrwrh(Sc, Sa, t(adjAC))  # size: n.A * n.C
  return( (mt.A2C + t(mt.C2A)) / 2 ) # size: n.C * n.A
}
