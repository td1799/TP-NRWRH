# TP-NRWRH 10 cross validation

source("utils.R", encoding = "UTF-8")
source("crossValidation.R", encoding = "UTF-8")
source("nrwrh.R", encoding = "UTF-8")

# select dataset to run cross validation on.
PREDICT.path <- "datasets/PREDICT/"
CDataset.path <- "datasets/CDataset/"

PREDICT.result.path <- "result_crossValidation/PREDICT/"
CDataset.result.path <- "result_crossValidation/CDataset/"

# on PREDICT as an example
data.path <- PREDICT.path

drugsName <- as.vector(read.table(paste0(data.path, "DrugsName"), header = F)[[1]])
diseasesName <- as.vector(read.table(paste0(data.path, "DiseasesName"), header = F)[[1]])

drugSim <- as.matrix(read.csv(paste0(data.path, "DrugSimMat"), header = F, sep = " "))
dimnames(drugSim) <- list(drugsName, drugsName)
diseaseSim <- as.matrix(read.csv(paste0(data.path, "DiseaseSimMat"), header = F, sep = " "))
dimnames(diseaseSim) <- list(diseasesName, diseasesName)

adjDiDr <- as.matrix(read.csv(paste0(data.path, "DiDrAMat"), header = F, sep = " "))
dimnames(adjDiDr) <- list(diseasesName, drugsName)

adjDrDi <- t(adjDiDr)

#参数设置
k.fold <- 10
LIMIT_TOR <- 1e-10
lambda <- 0.8
eta <- 0.4
restart.ratio <- 0.3

# about Zhou Tao 2010. In our method, this parameter is always 1.0.
flags.lambda <- 1.0

result.path <- PREDICT.result.path
if (!file.exists(file.path(".", result.path))) {
  dir.create(file.path(".", result.path), recursive = T)
}

result <- crossValidation(drugSim, diseaseSim, adjDrDi, k.fold)
auc <- computeAUC(c(result$scores), c(result$labels))
write.table(result$scores, paste0(result.path, "TP-NRWRH_scores.txt"), row.names = F, col.names = F, sep = " ")
write.table(result$labels, paste0(result.path, "TP-NRWRH_labels.txt"), row.names = F, col.names = F, sep = " ")
print(paste0("two-pass (TP-NRWRH) auc: ", auc))

result <- crossValidation.onePass(drugSim, diseaseSim, adjDrDi, k.fold)
auc <- computeAUC(c(result$scores), c(result$labels))
write.table(result$scores, paste0(result.path, "scores_Dr2Di.txt"), row.names = F, col.names = F, sep = " ")
write.table(result$labels, paste0(result.path, "labels_Dr2Di.txt"), row.names = F, col.names = F, sep = " ")
print(paste0("Dr2Di auc ", auc))

result <- crossValidation.onePass(diseaseSim, drugSim, adjDiDr, k.fold)
auc <- computeAUC(c(result$scores), c(result$labels))
write.table(result$scores, paste0(result.path, "scores_Di2Dr.txt"), row.names = F, col.names = F, sep = " ")
write.table(result$labels, paste0(result.path, "labels_Di2Dr.txt"), row.names = F, col.names = F, sep = " ")
print(paste0("Di2Dr auc ", auc))
