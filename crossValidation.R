source("TP-NRWRH.r", encoding = "UTF-8")

crossValidation <- function(Sa, Sc, adjAC, k.fold=10) {
  n.A <- nrow(Sa); n.C <- nrow(Sc)
  
  folds <- 1:k.fold
  scores <- matrix(0, nrow = n.A, ncol = n.C)
  masker <- matrix(sample(folds, n.A * n.C, T), nrow = n.A, ncol = n.C)
  for (test.fold in folds) {
    tmp.adjAC <- adjAC
    tmp.adjAC[which(masker == test.fold)] <- 0
    
    mt <- TP_NRWRH(Sa, Sc, tmp.adjAC)
    
    scores[which(masker == test.fold)] <- mt[which(masker == test.fold)]
  }
  return(list(
    scores    =   scores,
    labels    =   adjAC
  ))
}

crossValidation.onePass <- function(Sa, Sc, adjAC, k.fold=10) {
  n.A <- nrow(Sa); n.C <- nrow(Sc)
  
  folds <- 1:k.fold
  scores <- matrix(0, nrow=n.A, ncol=n.C)
  masker <- matrix(sample(folds, n.A * n.C, T), nrow=n.A, ncol=n.C)
  for (test.fold in folds) {
    tmp.adjAC <- adjAC
    tmp.adjAC[which(masker == test.fold)] <- 0
    
    mt <- TP_NRWRH.onePass(Sa, Sc, tmp.adjAC)
    
    mt <- t(mt)
    scores[which(masker == test.fold)] <- mt[which(masker == test.fold)]
  }
  return(list(
    scores    =   scores,
    labels    =   adjAC
  ))
}
