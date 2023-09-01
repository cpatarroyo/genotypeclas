library(caret)
library(poppr)
library(nnet)
library(adabag)

setwd("~/Uniandes/Rfiles/genotypeclas")

poptrain <- read.genalex("Training_DB.csv", ploidy = 3)

modTrain <- function(population) {
  
  #Prepare the population table for model training
  population <- clonecorrect(population, strata = NA)
  population <- missingno(population, cutoff = 0, type = "geno")
  poptrain <- population@tab
  poplabel <- population@pop
  
  #Remove uninformative "alleles" 
  uninf <- grep("\\.0$",colnames(poptrain))
  poptrain <- poptrain[,-uninf]
  
  #Train the neural network model
  set.seed(999)
  #searchspace <- expand.grid(size =1:8, decay = seq(0,5, by =0.5))
  #trmodel <- train(poptrain,poplabel, method = 'nnet',tuneLength = 1, tuneGrid = searchspace)
  trmodel <- caret::train(poptrain,poplabel, method = 'AdaBag', maxdepth = 30)
  
  #Save and return the trained model
  saveRDS(trmodel, file = "clasificator/Trained_model.rds")
  print(c((mean(trmodel$resample$Accuracy)),mean(trmodel$resample$Kappa)))
  #return(trmodel)
}

modTrain(poptrain)
