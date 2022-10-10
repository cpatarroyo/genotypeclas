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
  #Load the population to be classified
  newpop <- read.genalex(newdata,ploidy =3)
  
  #Load the trained model 
  trainedModel <- readRDS("clasificator/Trained_model.rds")
  
  #Load the allele table 
  newtable <- newpop@tab
  
  #Remove uninformative alleles
  uninf <- grep("\\.0$",colnames(newtable))
  newtable <- newtable[,-uninf]
  
  #Equate the predictors in the model to the variables in the new dataset
  missingAl <- which(!trainedModel$finalModel$xNames %in% colnames(newtable))
  tempMat <- matrix(data = 0, nrow = length(newtable[,1]), ncol = length(missingAl), dimnames = list(rownames(newtable), trainedModel$finalModel$xNames[missingAl]))
  newtable <- cbind(newtable,tempMat)

  #Generate the prediction and arrange in a results table
  resTable <- data.frame(Genotype = predict(trainedModel,newtable), Probability = unlist(apply(predict(trainedModel,newtable,type = "prob"),MARGIN = 1,FUN = max)))
  
  #Return the resulting table
  return(resTable)
  
}

test <- newPredict("Base de datos global_Camilo-CHN.csv")
