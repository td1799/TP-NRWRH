# compute top Rank of cutoffs, on PREDICT as an example

source("utils.R", encoding = "UTF-8")
source("nrwrh.R", encoding = "UTF-8")
source("TP-NRWRH.r", encoding = "UTF-8")

PREDICT.path <- "datasets/PREDICT/"
CDataset.path <- "datasets/CDataset/"

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

# set parameters
k.fold <- 10
LIMIT_TOR <- 1e-10 # convergence

lambda <- 0.8
eta <- 0.4
restart.ratio <- 0.3

# lambda in zhouTao 2012. always set as 1.0 in our work
flags.lambda <- 1.0

topRank <- computeTopRank(diseaseSim, drugSim, adjDiDr, k.fold)

print(topRank)
