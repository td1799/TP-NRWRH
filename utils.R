library(parallel)
library(ROCR)

parallel.matrix.multiply <- function(A, B) {
  if (ncol(A) != nrow(B)) {
    stop("Matrix A and B, not confirm.")
  }
  nr <- detectCores() - 2
  cl <- makeCluster(rep("localhost", nr))
  t0 <- proc.time()
  idx <- splitIndices(nrow(A), length(cl))
  Alist <- lapply(idx, function(ii) {
    return(A[ii,,drop=F])
  })
  ans <- clusterApply(cl, Alist,
    function(aa, BB) {
      return(aa %*% BB)  
    }, B
  )
  t1 <- proc.time()
  stopCluster(cl)
  return(do.call(rbind, ans))
}

laplacian_norm <- function(mtx) {
  D <- rowSums(mtx)
  C <- colSums(mtx)
  if (!isTRUE(all.equal(D, C, tolerance = 0.0001))) {
    print("laplacian norm error.")
  }
  D <- diag(1.0 / D)
  dimnames(D) <- list(rownames(mtx), rownames(mtx))
  D[is.infinite(D) | is.na(D)] <- 0
  D <- sqrt(D)
  return(D %*% mtx %*% D)
}

# 资源转移量，计算邻接矩阵A的行元素间的w[i][j]
# w[i][j] -> 如果已知一个列元素与row[i]相关，则它与row[j]相关的似然
# 返回 -> 一个W[n][n]的矩阵
# lambda=1时，完全考虑出边的平均；lambda=0时完全考虑入边的平均；
transferWeight <- function(n, m, A, lambda = 1.0) {
  if (nrow(A) != n | ncol(A) != m) {
    stop(paste0("传入的矩阵A应该为n*m的矩阵.", " n: ", n, " m: ", m ))
  }
  ky <- colSums(A)
  ky <- diag(1/ky) # 矩阵乘法用列向量相乘来理解, m*m
  kx <- rowSums(A) # length: n
  ky[is.infinite(ky) | is.na(ky)]	 <- 0 # Bugfix: 1/0=Infinite replaced with 0
  
  fx = 1 / (matrix(kx, nrow=n, ncol=n, byrow=FALSE)^(1-lambda) * 
              matrix(kx, nrow=n, ncol=n, byrow=TRUE)^(lambda)) #对应元素逐一相乘
  fx[is.infinite(fx) | is.na(fx)] <- 0
  
  W <- fx * ( (A %*% ky) %*% t(A) )
  return (W)
}

# computed new similarity matrix among X.
# return:
#   Sx    :      Sx(i, j): the number of known Ys shared by Xi and Xj.
computeSharedSn <- function(adjXY) {
  return(adjXY %*% t(adjXY))
}

L1_norm <- function(v) {
  return( sum(abs(v)) ) 
}

computeAUC <- function(scores, labels) {
  pred <- prediction(scores, labels)
  perf <- performance(pred, "tpr", "fpr")
  auc <- performance(pred, "auc")
  auc <- auc@y.values[[1]]
  return(auc)
}

sortWithinCol <- function(scoreMtx, labelMtx) {
  result <- matrix(0, nrow = nrow(scoreMtx), ncol = ncol(labelMtx))
  for (i in 1:ncol(scoreMtx)) {
    result[,i] <- labelMtx[,i][order(scoreMtx[,i],decreasing = T)]
  }
  return(result)
}

computeTopRank <- function(Sa, Sc, adjAC, k.fold=10) {
  A <- rownames(Sa); C <- colnames(Sc)
  n.A <- length(A); n.C <- length(C)
  
  result <- matrix(0, nrow=1, ncol = 5)
  cutoff <- c(1, 10, 20, 50, 100)
  
  folds <- 1:k.fold
  pairs <- which(adjAC == 1, arr.ind = T)
  n.pairs <- nrow(pairs)
  id <- sample(folds, n.pairs, replace = T)
  for (test.fold in folds) {
    test.set <- subset(pairs, id %in% test.fold)
    train.set <- subset(pairs, id %in% folds[-test.fold])
    tmp.adjAC <- adjAC
    if ( !all(tmp.adjAC[test.set] == 1) ) {
      stop("test tmp.adjAC error.")
    }
    tmp.adjAC[test.set] <- 0
    tmp.Sa <- Sa
    tmp.Sc <- Sc
    
    mt <- TP_NRWRH(tmp.Sa, tmp.Sc, tmp.adjAC)
    
    scoreMtx <- mt
    scoreMtx[train.set] <- 0
    
    labelMtx <- matrix(0, nrow = nrow(tmp.adjAC), ncol = ncol(tmp.adjAC))
    labelMtx[test.set] <- 1
    
    panel <- sortWithinCol(scoreMtx = scoreMtx, labelMtx = labelMtx)
    
    for (i in 1:length(cutoff)) {
      result[1,i] <- result[1,i] + sum(panel[1:cutoff[i], ])
    }
  }
  return(result)
}

computeTopRank_indepTest <- function(scoreMtx, labelMtx) {
  result <- matrix(0, nrow=1, ncol = 5)
  cutoff <- c(1, 10, 20, 50, 100)
  panel <- sortWithinCol(scoreMtx = scoreMtx, labelMtx = labelMtx)
  for (i in 1:length(cutoff)) {
    result[1,i] <- result[1,i] + sum(panel[1:cutoff[i], ])
  }
  return(result)
}
