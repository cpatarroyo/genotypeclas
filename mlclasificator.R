#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
args <- c("nn", "lolz")

library(RSNNS)
library(bnclassify)
library(klaR)
library(poppr)
library(nnet)
library(arm)
library(adabag)
library(plyr)
library(bst)
library(plyr)
library(obliqueRF)
library(party)
library(caret)

setwd("~/Uniandes/Rfiles/genotypeclas")

population <- read.genalex("Balanced_DS.csv", ploidy = 3)

modTrain <- function(population, method = NULL) {
  
  #Prepare the population table for model training
  population <- clonecorrect(population, strata = NA)
  population <- missingno(population, cutoff = 0, type = "geno")
  poptrain <- population@tab
  poplabel <- population@pop
  
  #Remove uninformative "alleles" 
  uninf <- grep("\\.0$",colnames(poptrain))
  poptrain <- poptrain[,-uninf]
  
  set.seed(999)
  
  #Training method selection
  if(method == "nn" || is.null(method)) {
    #Train the neural network model
    searchspace <- expand.grid(size =1:8, decay = seq(0,5, by =0.5))
    trmodel <- caret::train(poptrain,poplabel, method = 'nnet',tuneLength = 1, tuneGrid = searchspace)
  }
  else if(method == "nb") {
    #Train the naive bayes clasificator model
    trmodel <- caret::train(poptrain,poplabel, method = 'nb',tuneLength = 1)
  }
  else if(method == "bayesglm") {
    #Train Model by a bayes generalized linear model
    trmodel <- caret::train(poptrain,poplabel, method = 'bayesglm')
  }
  else if(method == "adabag") {
    #Train model by Ada boost
    trmodel <- caret::train(poptrain,poplabel, method = 'AdaBag', maxdepth = 30)
  } 
  else if(method == "bsttree") {
    #Train model by Ada boost classification tree
    trmodel <- caret::train(poptrain,poplabel, method = 'bstTree')
  }
  else if(method == "oblique") {
    trmodel <- caret::train(poptrain,poplabel, method = 'ORFpls')
  }
  else if(method == "cforest") {
    trmodel <- caret::train(poptrain,poplabel, method = 'cforest')
  }
  else if(method == "mlnn") {
    searchspace <- expand.grid(layer1=1:5,layer2=0:5,layer3=0:5)
    trmodel <- caret::train(poptrain,poplabel, method = 'mlpML',tuneGrid = searchspace)
  }
  else {
    stop("You must enter a valid training method")
  }
  
  #Save and return the trained model
  #saveRDS(trmodel, file = "clasificator/Trained_model.rds")
  print(c((mean(trmodel$resample$Accuracy)),mean(trmodel$resample$Kappa)))
  return(trmodel)
}

newPredict <- function(newdata, model = NULL) {
  #Load the population to be classified
  #newpop <- read.genalex(newdata,ploidy =3)
  newpop <- newdata
  
  #Load the trained model 
  if(is.null(model)) {
    trainedModel <- readRDS("clasificator/Trained_model.rds")
  } else {
    trainedModel <- model
  }
  
  #Load the allele table 
  newtable <- newpop@tab
  
  #Remove uninformative alleles
  uninf <- grep("\\.0$",colnames(newtable))
  newtable <- newtable[,-uninf]
  
  if(class(trainedModel$finalModel)[1] != "RandomForest") {
    #Equate the predictors in the model to the variables in the new dataset
    missingAl <- which(!trainedModel$finalModel$xNames %in% colnames(newtable))
    tempMat <- matrix(data = 0, nrow = length(newtable[,1]), ncol = length(missingAl), dimnames = list(rownames(newtable), trainedModel$finalModel$xNames[missingAl]))
    newtable <- cbind(newtable,tempMat)
  }
  
  #Generate the prediction and arrange in a results table
  resTable <- data.frame(Name = rownames(newtable),
                         PrLineage = predict(trainedModel,newtable),
                         #Probability = unlist(apply(predict(trainedModel,newtable,type = "prob"),MARGIN = 1,FUN = max)),
                         Lineage = newpop@pop)
  
  #Return the resulting table
  return(resTable)
  
}

count <- 0
restab <- data.frame(Accuracy = double(), Kappa = double(), TestPred = double())

#Niter tests for the accuracy of the ML algorithm prediction
while(count < 3) {
  #Define the data partition for training and testing
  trainlen <- round(summary(population)$n*0.2,digits = 0)
  trainindex <- sort(sample(1:summary(population)$n,trainlen,replace = FALSE))
  training <- population[-trainindex]
  test <- population[trainindex]
  test <- missingno(test, cutoff = 0, type = "geno")
  
  #Train the ML model
  trainedModel <- modTrain(training, method = args[1])
  
  #Make the prediction
  prediction <- newPredict(test,model = trainedModel)
  restab <- rbind(restab,c(trainedModel$results$Accuracy[tolerance(trainedModel$results,metric = "Accuracy",maximize = TRUE)],trainedModel$results$Kappa[tolerance(trainedModel$results,metric = "Kappa",maximize = TRUE)],sum(as.character(prediction$PrLineage) == as.character(prediction$Lineage),na.rm = T)/(dim(prediction)[1])))
  count <- count+1
}

colnames(restab)<-c("Accuracy","Kappa","TestAc")

write.csv(restab,file = paste(args[2],".csv",sep = ""))