install.packages("tidyverse")
install.packages("dplyr")
install.packages("tidyr")
install.packages("devtools")
install.packages("DSR")
install.packages("gridExtra")
install.packages('caTools')
library(caTools)
library(rpart.plot)
library(rattle)
library(rpart)
library(heuristica)
library(caret)
library(rpart.plot)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DSR)
library(party)

cancerData<-read.csv(file.choose())
head(cancerData)
summary(cancerData)

#Insert NA's data
cancerData$Bare.nuclei <- ifelse(is.na(cancerData$Bare.nuclei), 
                                 round(ave(cancerData$Bare.nuclei, 
                                           FUN = function(x) 
                                             mean(x, na.rm = TRUE ))),
                                 cancerData$Bare.nuclei)

head(cancerData)

#Delete Id (Make the data easier)
cancerData <- select(cancerData, -Id) 
head(cancerData)

#Convert to factor
cancerData$Cl.thickness = factor(cancerData$Cl.thickness)
cancerData$Cell.size = factor(cancerData$Cell.size)
cancerData$Cell.shape = factor(cancerData$Cell.shape)
cancerData$Marg.adhesion = factor(cancerData$Marg.adhesion)
cancerData$Epith.c.size = factor(cancerData$Epith.c.size)
cancerData$Bare.nuclei = factor(cancerData$Bare.nuclei)
cancerData$Bl.cromatin = factor(cancerData$Bl.cromatin)
cancerData$Normal.nucleoli = factor(cancerData$Normal.nucleoli)
cancerData$Mitoses = factor(cancerData$Mitoses)
cancerData$Class = factor(cancerData$Class)

#check the class (make sure the output is 'factor')
class(cancerData$Cl.thickness)
class(cancerData$Cell.size)
class(cancerData$Cell.shape)
class(cancerData$Marg.adhesion)
class(cancerData$Epith.c.size)
class(cancerData$Bare.nuclei)
class(cancerData$Bl.cromatin)
class(cancerData$Normal.nucleoli)
class(cancerData$Mitoses)
class(cancerData$Class)

#Random
set.seed(seed<-489)

#Split (Data Training & Testing)
split = sample.split(cancerData$Class, SplitRatio = 0.8)
cancerTrainingData = subset(cancerData, split == TRUE)
cancerTestData = subset(cancerData, split == FALSE)

#check dimension data
dim(cancerTrainingData)
dim(cancerTestData)

prop.table(table(cancerTrainingData$Class))
prop.table(table(cancerTestData$Class))

#MODEL: rpart
cancerModelR <- rpart(Class~., data = cancerTrainingData, 
                      method = 'class')
summary(cancerModelR)
cancerModelR$variable.importance
barplot(cancerModelR$variable.importance)
rpart.plot(cancerModelR)

#MODEL: ctree
cancerModelC <- ctree(Class ~ ., data = cancerTrainingData)
plot(cancerModelC)

#Predict: Use ctree
predictCancerModelC = predict(cancerModelC, 
                              newdata = cancerTestData)
table(predictCancerModelC, cancerTestData$Class)
confusionMatrix(data=predictCancerModelC,reference=cancerTestData$Class)

#Predict: Use rpart
predictCancerModelR = predict(cancerModelR, 
                              newdata = cancerTestData, 
                              type = "class")
summary(predictCancerModelR)
table(predictCancerModelR, cancerTestData$Class)
confusionMatrix(data=predictCancerModelR,reference=cancerTestData$Class)