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
library("rstudioapi")

# Data loading
name <- "Japonica"
Dir <- gsub("/RF_GY.R", "", getSourceEditorContext()$path)
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
length(rownames(Markers)) == length(unique(rownames(Markers)))#320 l?neas
unique(Pheno$Line) %>% length()
Geno = tcrossprod(as.matrix(Markers[, 2:16384]))/dim(Markers)[2]
dim(Geno)

# Data preparation of Env, G & GE
Line <- model.matrix(~0 + Line, data = Pheno)
Env <- model.matrix(~0 + Env, data = Pheno)
Geno <- t(chol(Geno)) 
LineG <- Line %*% Geno
LinexGenoxEnv <- model.matrix(~ 0 + LineG:Env)

# Predictor and Response Variables
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

# Model training and predictions of the ith partition
for (i in seq(folds)) {
  cat("*** Fold:", i, " ***\n")
  # Identify the training and testing sets
  fold <- folds[[i]]
  
  # Identify the training and testing sets
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training]
  y_testing <- y[fold$testing]
  
  # Model training
  model <- random_forest(
    x = X_training,
    y = y_training,
    
    # Specify the hyperparameters ranges
    trees_number = list(min = 50, max = 500),
    node_size = list(min = 5, max =15),
    
    tune_type = "Bayesian_optimization",
    tune_bayes_samples_number = 5,
    tune_bayes_iterations_number = 5,
    
    #The iterations wont be shown
    verbose = FALSE
  )
  
  # Testing Predictions
  predictions <- predict(model, X_testing)
  
  # Predictions for the Variable 1: GY
  FoldPredictionsV1 <- data.frame(
    Fold = i,
    Line = Pheno$Line[fold$testing],
    Env = Pheno$Env[fold$testing],
    Observed = y_testing,
    Predicted = predictions$V1$predicted)
  PredictionsV1 <- rbind(PredictionsV1, FoldPredictionsV1)
  
  # Hyperparams
  HyperparamsFold <- model$hyperparams_grid %>%
    mutate(Fold = i)
  Hyperparams <- rbind(Hyperparams, HyperparamsFold)
  
  # Best hyperparams
  print(model$best_hyperparams)
}

#######################
###### Summaries ######
#######################
res <- "RF"
mkdir(paste(name, res, sep = ""))

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

# Summaries by Environment
SUM_tst_env_V1=summariesV1$env
data.table::fwrite(
  SUM_tst_env_V1,
  file = paste(name, res, "/Summary_GY_testing_E+G+GxE.csv", sep = ""),
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)

#Hyperparams
data.table::fwrite(
  Hyperparams,
  file = paste(name, res, "/Hyperparams.csv", sep = ""),
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)