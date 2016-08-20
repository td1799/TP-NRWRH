# do independent test
# train on PREDICT, and test on another dataset from MBiRW

source("utils.R", encoding = 'utf-8')
source("nrwrh.R", encoding = 'utf-8')
source("TP-NRWRH.r", encoding = 'utf-8')

PREDICT.path <- "datasets/PREDICT/"

train.adjDiDr <- as.matrix(read.table(paste0(PREDICT.path, "DiDrAMat"), header = F, sep = ' '))
train.drugsName <- as.vector(read.csv(paste0(PREDICT.path, "DrugsName"), header = F)[[1]])
train.diseasesName <- as.vector(read.csv(paste0(PREDICT.path, "DiseasesName"), header = F)[[1]])
dimnames(train.adjDiDr) <- list(train.diseasesName, train.drugsName)
train.drugSim <- as.matrix(read.table(paste0(PREDICT.path, "DrugSimMat"), header = F, sep = " "))
train.diseaseSim <- as.matrix(read.table(paste0(PREDICT.path, "DiseaseSimMat"), header = F, sep = " "))
dimnames(train.drugSim) <- list(train.drugsName, train.drugsName)
dimnames(train.diseaseSim) <- list(train.diseasesName, train.diseasesName)

# parameters
LIMIT_TOR <- 1e-10
lambda <- 0.7
eta <- 0.5
restart.ratio <- 0.7
flags.lambda <- 1.0

result.path <- "result_indepTest/"
if (!file.exists(file.path("./", result.path))) {
  dir.create(file.path("./", result.path), recursive = T)
}

predictions <- TP_NRWRH(train.diseaseSim, train.drugSim, train.adjDiDr)
write.table(predictions, paste0(result.path, "TP-NRWRH_indepTest_prediction.txt"), row.names = F, col.names = F, sep = ' ')

# compute top rank

TP_NRWRH.predictions <- "TP-NRWRH_indepTest_prediction.txt"

predictions <- as.matrix(read.table(paste0(result.path, TP_NRWRH.predictions), header = F, sep = ' '))

rownames(predictions) <- as.vector(read.csv(paste0(PREDICT.path, "DiseasesName"), header = F)[[1]])
colnames(predictions) <- as.vector(read.csv(paste0(PREDICT.path, "DrugsName"), header = F)[[1]])

train <- as.matrix(read.table(paste0(PREDICT.path, "DiDrAMat")))
rownames(train) <- as.vector(read.csv(paste0(PREDICT.path, "DiseasesName"), header = F)[[1]])
colnames(train) <- as.vector(read.csv(paste0(PREDICT.path, "DrugsName"), header = F)[[1]])

indepTest.path = "datasets/indepTest/"
test <- as.matrix(read.table(paste0(indepTest.path, "DiDrAMat")))
rownames(test) <- as.vector(read.csv(paste0(indepTest.path, "DiseasesName"), header = F)[[1]])
colnames(test) <- as.vector(read.csv(paste0(indepTest.path, "DrugsName"), header = F)[[1]])

predictions <- predictions[rownames(test), colnames(test)]
train <- train[rownames(test), colnames(test)]

predictions[train == 1] <- 0
labels <- test

# cutoff: 1 10 20 50 100
resultTop <- computeTopRank_indepTest(scoreMtx = predictions, labelMtx = labels)
print(resultTop)

scores <- c(predictions[train != 1])
labels <- c(test[train != 1])
auc <- computeAUC(scores, labels)
print(auc)
