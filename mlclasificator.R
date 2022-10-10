library(caret)
#library(bnclassify)
#library(klaR)
library(poppr)
library(nnet)

setwd("C:/Users/capat/Documents/Uniandes/Tesis/Tangent/genotypeclas")

poptrain <- read.genalex("Base de datos ColPeru.csv", ploidy = 3)

modTrain <- function(population) {
  
  #Prepare the population table for model training
  #population <- read.genalex(population)
  population <- clonecorrect(population, strata = NA)
  population <- missingno(population, cutoff = 0, type = "geno")
  poptrain <- population@tab
  poplabel <- population@pop
  
  #Remove uninformative "alleles" 
  uninf <- grep("\\.0$",colnames(poptrain))
  poptrain <- poptrain[,-uninf]
  
  #Train the neural network model
  set.seed(999)
  searchspace <- expand.grid(size =1:8, decay = seq(0,5, by =0.5))
  trmodel <- train(poptrain,poplabel, method = 'nnet',tuneLength = 1, tuneGrid = searchspace)
  #Train the naive bayes clasificator model
  #trmodel <- train(poptrain,poplabel, method = 'nb',tuneLength = 1)
  
  #Save and return the trained model
  saveRDS(trmodel, file = "Trained_model.rds")
  print(c((mean(trmodel$resample$Accuracy)),mean(trmodel$resample$Kappa)))
  return(trmodel)
}

newPredict <- function(newdata) {
  poptest <- read.genalex("population.csv",ploidy =3)
}

testmod <- modtrain(population = population)

confusionMatrix.train(trainedmod)

trainedmod$finalModel$xNames
