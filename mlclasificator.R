library(caret)
#library(bnclassify)
#library(klaR)
library(poppr)
library(nnet)

setwd("C:/Users/capat/Documents/Uniandes/Tesis/Tangent/genotypeclas")

modtrain <- function(population) {
  
  #Prepare the population table for model training
  #population <- read.genalex(population)
  population <- clonecorrect(population, strata = NA)
  population <- missingno(population, cutoff = 0, type = "geno")
  poptrain <- population@tab
  poplabel <- population@pop
  
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
  
}

testmod <- modtrain(population = population)

confusionMatrix.train(testmod)

testmod$finalModel$xNames
