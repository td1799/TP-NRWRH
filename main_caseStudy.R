# drug-centered & disease-centered (one-pass & two-pass) case study

source("utils.R", encoding = "UTF-8")
source("nrwrh.R", encoding = "UTF-8")
source("TP-NRWRH.r", encoding = "UTF-8")

PREDICT.path <- "datasets/PREDICT/"

# on PREDICT
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
eta <- 0.3
restart.ratio <- 0.3
flags.lambda <- 1.0

result.path <- "result_caseStudy/"
if (!file.exists(file.path(".", result.path))) {
  dir.create(file.path(".", result.path), recursive = T)
}

mt <- TP_NRWRH(diseaseSim, drugSim, adjDiDr)
print(paste("Size of adjDiDr:", nrow(adjDiDr), ncol(adjDiDr)))
print(paste("Size of mt:", nrow(mt), ncol(mt)))
if (nrow(mt) != nrow(adjDiDr) || ncol(mt) != ncol(adjDiDr)) {
  stop("avg drug error")
}
if (!all(rownames(mt) == rownames(adjDiDr)) || !all(colnames(mt) == colnames(adjDiDr))) {
  stop("avg drug name error")
}
mt[adjDiDr == 1] <- 0
result <- list()
for (j in 1:ncol(mt)) {
  result[[ colnames(mt)[j] ]] <- sort(mt[,j], decreasing = T)
}
saveRDS(result, paste0(result.path, 'caseStudy_avg_aboutDrug.RDS'))


mt <- TP_NRWRH(drugSim, diseaseSim, adjDrDi)
print(paste("Size of adjDrDi: ", nrow(adjDrDi), ncol(adjDrDi)))
print(paste("Size of mt:", nrow(mt), ncol(mt)))
if (nrow(mt) != nrow(adjDrDi) || ncol(mt) != ncol(adjDrDi)) {
  stop("avg disease error")
}
if (!all(rownames(mt) == rownames(adjDrDi)) || !all(colnames(mt) == colnames(adjDrDi))) {
  stop("avg disease name error.")
}
mt[adjDrDi == 1] <- 0 # drug - disease
result <- list()
for (j in 1:ncol(mt)) {
  result[[ colnames(mt)[j] ]] <- sort(mt[,j], decreasing = T)
}
saveRDS(result, paste0(result.path, 'caseStudy_avg_aboutDisease.RDS'))

mt <- TP_NRWRH.onePass(diseaseSim, drugSim, adjDiDr)
print(paste("Size of adjDiDr:", nrow(adjDiDr), ncol(adjDiDr)))
print(paste("Size of mt:", nrow(mt), ncol(mt)))
if (nrow(mt) != ncol(adjDiDr) || ncol(mt) != nrow(adjDiDr)) {
  stop("oneWay disease error")
}
if (!all(rownames(mt) == colnames(adjDiDr)) || !all(colnames(mt) == rownames(adjDiDr))) {
  stop("oneWay disease names error.")
}
mt[t(adjDiDr) == 1] <- 0
result <- list()
for (j in 1:ncol(mt)) {
  result[[ colnames(mt)[j] ]] <- sort(mt[,j], decreasing = T)
}
saveRDS(result, paste0(result.path, 'caseStudy_onePass_aboutDisease_fromDisease2Drug.RDS'))


mt <- TP_NRWRH.onePass(drugSim, diseaseSim, adjDrDi)
print(paste("Size adjDrDi:", nrow(adjDrDi), ncol(adjDrDi)))
print(paste("Size mt:", nrow(mt), ncol(mt)))
if (nrow(mt) != ncol(adjDrDi) || ncol(mt) != nrow(adjDrDi)) {
  stop("oneWay drug error.")
}
if (!all(rownames(mt) == colnames(adjDrDi)) || !all(colnames(mt) == rownames(adjDrDi))) {
  stop("oneWay drug name error.")
}
mt[t(adjDrDi) == 1] <- 0
result <- list()
for (j in 1:ncol(mt)) {
  result[[ colnames(mt)[j] ]] <- sort(mt[,j], decreasing = T)
}
saveRDS(result, paste0(result.path, 'caseStudy_onePass_aboutDrug_fromDrug2Disease.RDS'))


mt <- TP_NRWRH.onePass(drugSim, diseaseSim, adjDrDi)
print(paste("Size adjDrDi:", nrow(adjDrDi), ncol(adjDrDi)))
print(paste("Size mt:", nrow(mt), ncol(mt)))
if (nrow(mt) != ncol(adjDrDi) || ncol(mt) != nrow(adjDrDi)) {
  stop("oneWay drug error.")
}
if (!all(rownames(mt) == colnames(adjDrDi)) || !all(colnames(mt) == rownames(adjDrDi))) {
  stop("oneWay drug name error.")
}
mt[t(adjDrDi) == 1] <- 0
mt <- t(mt)
result <- list()
for (j in 1:ncol(mt)) {
  result[[ colnames(mt)[j] ]] <- sort(mt[,j], decreasing = T)
}
saveRDS(result, paste0(result.path, 'caseStudy_onePass_aboutDisease_fromDrug2Disease.RDS'))

result <- readRDS(paste0(result.path, "caseStudy_avg_aboutDisease.RDS"))
result$D605055[1:20]

