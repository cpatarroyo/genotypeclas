library(caret)
library(bnclassify)
library(klaR)
library(poppr)

population <- read.genalex("Base de datos ColPeru.csv",ploidy = 3)


#bruvo.msn(population, vertex.label = NA)
population <- clonecorrect(population, strata = NA)

population <- missingno(population, cutoff = 0, type = "geno")

str(population)

population@tab
population@pop

restest <- data.frame(Real = population@pop, Pred = NA)
rownames(restest) <- rownames(population@tab)

#head(restest)

training <- sort(sample(1:362,floor(0.7*length(population@pop))))

poptrain <- population@tab[training,]
poplabel <- as.character(population@pop[training])

poptest <- population@tab[-training,]
restest <- restest[-training,]

set.seed(999)
model <- train(poptrain, poplabel, method = 'nb',tuneLength = 10)
warnings()
model

restest$Pred <- predict(model,newdata = poptest)
confusionMatrix.train(model)

