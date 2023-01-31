library(caret)
#library(bnclassify)
library(klaR)
library(poppr)
library(nnet)
#library(arm)
#library(adabag)
#library(plyr)

setwd("~/Uniandes/Rfiles/genotypeclas")

population <- read.genalex("Training_DB.csv", ploidy = 3)

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
  
  #Train Model by a bayes generalized linear model
  #trmodel <- train(poptrain,poplabel, method = 'bayesglm')
  
  #Train model by Ada boost
  #trmodel <- train(poptrain,poplabel, method = 'AdaBag', maxdepth = 30)
  
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
  
  #Equate the predictors in the model to the variables in the new dataset
  missingAl <- which(!trainedModel$finalModel$xNames %in% colnames(newtable))
  tempMat <- matrix(data = 0, nrow = length(newtable[,1]), ncol = length(missingAl), dimnames = list(rownames(newtable), trainedModel$finalModel$xNames[missingAl]))
  newtable <- cbind(newtable,tempMat)

  #Generate the prediction and arrange in a results table
  resTable <- data.frame(rownames(newtable),
                         PrLineage = predict(trainedModel,newtable),
                         Probability = unlist(apply(predict(trainedModel,newtable,type = "prob"),MARGIN = 1,FUN = max)),
                         Lineage = newpop@pop)
  
  #Return the resulting table
  return(resTable)
  
}

count <- 0
restab <- data.frame(Accuracy = double(), Kappa = double(), TestPred = double())

#Niter tests for the accuracy of the ML algorithm prediction
while(count < 20) {
  #Define the data partition for training and testing
  trainlen <- round(summary(population)$n*0.8,digits = 0)
  trainindex <- sort(sample(1:summary(population)$n,trainlen,replace = FALSE))
  training <- population[trainindex]
  test <- population[-trainindex]
  
  #Train the ML model
  trainedModel <- modTrain(training)
  indexRes <- which(trainedModel$results$size==trainedModel$bestTune$size)[which(which(trainedModel$results$size==trainedModel$bestTune$size) %in% which(trainedModel$results$decay == trainedModel$bestTune$decay))]
  
  #Make the prediction
  prediction <- newPredict(test)
  restab <- rbind(restab,c(trainedModel$results$Accuracy[indexRes],trainedModel$results$Kappa[indexRes],sum(as.character(prediction$PrLineage) == as.character(prediction$Lineage),na.rm = T)/(dim(prediction)[1])))
  
}

write.csv(restab,file = "Output.csv")