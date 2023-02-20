#Process to select a balanced learning set

#Selecting the lineages with higher representation
test1eq <- popsub(population,sublist = c("EC1","EU2A1","EU13"))

repvec <- c(sample(1:1081,154,replace = FALSE),sample(which(pop(test1eq) == "EU2A1"),154, replace=FALSE),which(pop(test1eq) == "EU13"))

test1eq <- test1eq[repvec]

genind2genalex(test1eq,filename = "Balanced_DS.csv", overwrite = TRUE)

#Function to prepare the population Data for phytophtora ID

prepare <- function(route, index=NULL) {
  tempTab <- read.csv(route, header=T, skip = 2)
  mmsInd <- seq(3, length(tempTab[1,]), by=3)
  
  ssrTab <- data.frame(tempTab$Ind,"query")
  popTab <- data.frame(tempTab$Ind,tempTab$Pop)
  
  for(i in mmsInd) {
    ssrTab <- cbind(ssrTab, as.data.frame(paste(tempTab[,i],tempTab[,i+1],tempTab[,i+2],sep = "/")))
  }
  
  if(is.null(index)) {
    ssrname <- "Tabla_SSR.csv"
    popname <- "Real_pops.csv"
  }
  else {
    ssrname <- paste("Tabla_SSR",index,".csv",collapse = "_")
    popname <- paste("Real_pops",index,".csv",collapse = "_")
  }
  
  colnames(ssrTab) <- c("Id","Lineage",colnames(tempTab)[mmsInd])
  colnames(popTab) <- c("Id", "Lineage")
  
  write.csv(ssrTab,file = ssrname, row.names = FALSE)
  write.csv(popTab,file = popname, row.names = FALSE)
}
