rm(list = ls())
library(dplyr)
library(SKM)
library(pls)

name <- "Japonica"
Dir <- gsub("/predictions.R", "", getSourceEditorContext()$path)
setwd(Dir)

load(paste(name, ".RData", sep = ""), verbose = TRUE) #Cargar los datos
ls()

# Markers
dim(Markers)
colnames(Markers)

# Pheno
dim(Pheno)
colnames(Pheno)
unique(Pheno$Env)
length(unique(Pheno$Line))

# ¿De cuántos híbridos contamos con su genotipo?
Pheno %>%
  group_by(Env) %>%
  select(Line) %>%
  count()

New_Pheno <- data.frame(Env = c(rep("2009", 320), 
                                rep("2010", 320), 
                                rep("2011", 320), 
                                rep("2012", 320),
                                rep("2013", 320)),
                        Line = rep(unique(Pheno$Line), 5),
                        GY = NA
           )

New_Pheno = Pheno[, 1:3]

# Agregamos la lineas que no tienen genotipo
for(e in unique(Pheno$Env)){
  #e = 2009
  New_Lines = which(!unique(Pheno$Line) %in% Pheno[Pheno$Env ==e, ]$Line)
  unique(Pheno$Line)[New_Lines]
  
  New_hybrids <- data.frame(Env = e, 
                            Line = unique(Pheno$Line)[New_Lines],
                            GY = NA)
  New_Pheno = rbind(New_Pheno, New_hybrids)
}
dim(New_Pheno)

#GRM
length(rownames(Markers)) == length(unique(rownames(Markers)))#320 l?neas
unique(Pheno$Line) %>% length()
Geno = tcrossprod(as.matrix(Markers[, 2:16384]))/dim(Markers)[2]
dim(Geno)

# Data preparation of Env, G & GE
Line <- model.matrix(~0 + Line, data = New_Pheno)
Env <- model.matrix(~0 + Env, data = New_Pheno)
Geno <- t(chol(Geno)) 
LineG <- Line %*% Geno
LinexGenoxEnv <- model.matrix(~ 0 + LineG:Env)

# Predictor and Response Variables
X <- cbind(Env, LineG, LinexGenoxEnv)
#y <- cbind(Pheno$GY, Pheno$PHR, Pheno$GC, Pheno$PH)
y <- New_Pheno$GY

# Conjunto de datos de entrenamiento
pos_newData= which(is.na(y))

# Identify the training and testing sets
X_training <- X[-pos_newData, ]
X_new <- X[pos_newData, ]
y_training <- y[-pos_newData]

# Model training
model <- random_forest(
  x = X_training,
  y = y_training,
  
  # Specify the hyperparameters ranges
  trees_number = list(min = 50, max = 500),
  node_size = list(min = 5, max =15),
  
  tune_type = "Bayesian_optimization",    
  #In this example the iterations wont be shown
  verbose = FALSE
)

# Testing Predictions
predictions <- predict(model, X_new)

head(predictions)

write.csv(predictions,file="predictions.csv")