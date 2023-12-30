#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Load the required libraries
library(caret)
library(poppr)
library(nnet)

#Define the prediction function
newPredict <- function(newdata) {
  #Load the population to be classified
  newpop <- read.genalex(newdata,ploidy =3)
  
  #Load the trained model 
  trainedModel <- readRDS("classifier/Trained_model.rds")
  
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
  resTable <- data.frame(Lineage = predict(trainedModel,newtable), Probability = unlist(apply(predict(trainedModel,newtable,type = "prob"),MARGIN = 1,FUN = max)))
  
  #Return the resulting table
  return(resTable)
  
}

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "Prediction.csv"
}

outputtable <- newPredict(args[1])

write.csv(outputtable, file = args[2])