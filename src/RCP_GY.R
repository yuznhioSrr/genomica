# Install SKM if it is necessary 
if(! "SKM" %in% rownames(installed.packages())){
  if (!require("devtools")) {
    install.packages("devtools")
  }
  devtools::install_github("gdlc/BGLR-R")
  devtools::install_github("rstudio/keras")
  devtools::install_github("rstudio/tensorflow")
  devtools::install_github("brandon-mosqueda/SKM")
}

# Install dplyr if it is necessary 
if(! "dplyr" %in% rownames(installed.packages())){
  install.packages("dplyr")
}

# Install pls if it is necessary 
if(! "pls" %in% rownames(installed.packages())){
  install.packages("pls")
}
# Install rstudioapi if it is necessary 
if(! "rstudioapi" %in% rownames(installed.packages())){
  install.packages("rstudioapi")
}

rm(list = ls())
library(dplyr)
library(SKM)
library(pls)
library("rstudioapi")

# Data loading
name <- "Japonica"
Dir <- gsub("/RCP_GY.R", "", getSourceEditorContext()$path)
setwd(Dir)
load(paste(name, ".RData", sep = ""), verbose = TRUE) #Cargar los datos
ls()

# Matrix of  Markers
dim(Markers)
colnames(Markers)

# Pheno: Matrix of Phenotypes 
dim(Pheno)
colnames(Pheno)

# GRM: Genomic Relationship Matrix
length(rownames(Markers)) == length(unique(rownames(Markers)))#320 líneas
unique(Pheno$Line) %>% length()
Geno = tcrossprod(as.matrix(Markers[, 2:16384]))/dim(Markers)[2]
dim(Geno)

# Data preparation of Env, G & GE
Line <- model.matrix(~0 + Line, data = Pheno)
Env <- model.matrix(~0 + Env, data = Pheno)
Geno <- t(chol(Geno)) 
LineG <- Line %*% Geno
LinexGenoxEnv <- model.matrix(~ 0 + LineG:Env)

# Predictor Variables 
X <- cbind(Env, LineG, LinexGenoxEnv)
y <- Pheno$GY # Response Variable 

# 5-fold Cross Validation indexes
set.seed(2022)
folds <- cv_kfold(records_number = nrow(X),
                   k = 5)

# PredictionsV1: data frame that will contain the variables:
##              (Number) Fold, Line, Env, 
##              (testing values) Observe & Predicted (values)
PredictionsV1 <- data.frame()
Hyperparams <- data.frame()

data <- cbind(y,X)

# Model training and predictions of the ith partition
for (i in seq(folds)) {
  cat("*** Fold:", i, " ***\n")
  # Identify the training and testing sets
  fold <- folds[[i]]
  
  # Model training
  model <- pcr(y ~ . , data = as.data.frame(data), validation = "CV", subset = fold$training)
  n = selectNcomp(model, method = "onesigma", plot = TRUE) #selection of components
  
  # Model taining with the One Standard Error Rule (method = "onesigma")
  ## optimal number of components selected
  model_f <- pcr(y ~ ., data = as.data.frame(data), validation = "CV", ncomp = n, subset = fold$training)
  X_tst <- data[fold$testing, -1]
 
  predictions <- predict(model_f, newdata = X_tst, ncomp = n)
 
  # Predictions for the Variable 1: GY
  FoldPredictionsV1 <- data.frame(
    Fold = i,
    Line = Pheno$Line[fold$testing],
    Env = Pheno$Env[fold$testing],
    Observed = data[fold$testing, 1],
    Predicted = predictions[,1,1])
  PredictionsV1 <- rbind(PredictionsV1, FoldPredictionsV1)
  
  # Best hyperparams
  print(n)
}

#######################
###### Summaries ######
#######################
res <- "PCR"
mkdir(paste(name, res, sep = "")) #New directory

# Prediction of Variable 1: GY
data.table::fwrite(
  PredictionsV1,
  file = paste(name, res, "/PredictionsGY.csv", sep = ""),
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

# Summaries
summariesV1 <- gs_summaries(PredictionsV1)

# Summaries by environment 
SUM_tst_env_V1=summariesV1$env
data.table::fwrite(
  SUM_tst_env_V1,
  file = paste(name, res, "/Summary_GY_testing_E+G+GxE.csv", sep = ""),
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)
