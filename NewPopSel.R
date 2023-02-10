#Process to select a balanced learning set

#Selecting the lineages with higher representation
test1eq <- popsub(population,sublist = c("EC1","EU2A1","EU13"))

repvec <- c(sample(1:1081,154,replace = FALSE),sample(which(pop(test1eq) == "EU2A1"),154, replace=FALSE),which(pop(test1eq) == "EU13"))

test1eq <- test1eq[repvec]

genind2genalex(test1eq,filename = "Balanced_DS.csv", overwrite = TRUE)
