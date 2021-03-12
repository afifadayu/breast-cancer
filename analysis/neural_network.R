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
library(tidyverse)
library(neuralnet)
library(GGally)
library(ISLR)
library("caret")
attach(Carseats)
library(tidyverse)
library("keras")
library(neuralnet)
library(Hmisc)
library(caret)

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

#Encoding Categorical Data
cancerData$Class <- factor(cancerData$Class, 
                           levels = c('malignant', 'benign'), 
                           labels = c(1, 0))

#Data Visualization

#1 ggpairs
ggpairs(cancerData, title = "Scatterplot Matrix of the Features of the Cancer Data")

#2 histogram tiap class
cancerData %>% gather(Atributes, Value, 2:10) %>%
  ggplot(aes (x=Value, fill=Atributes)) +
  geom_histogram(colour="grey")+
  facet_wrap(~Class)+ theme_bw()+
  labs(x="Values", y="Frequency",
       title="Cancer Data",
       subtitle="Histogram for each nerves") +
  theme(legend.title=element_blank(),
        legend.position="bottom")

#3 histogram per variabel
h_Cl.thickness<-cancerData %>% 
  group_by(Class,Cl.thickness)%>%
  ggplot(aes(x=Cl.thickness, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Cl.thickness", x="",y="")
h_Cell.size<-cancerData %>% 
  group_by(Class,Cell.size)%>%
  ggplot(aes(x=Cell.size, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Cell.size", x="",y="")
h_Cell.shape<-cancerData %>% 
  group_by(Class,Cell.shape)%>%
  ggplot(aes(x=Cell.shape, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Cell.shape", x="",y="")
h_Marg.adhesion<-cancerData %>% 
  group_by(Class,Marg.adhesion)%>%
  ggplot(aes(x=Marg.adhesion, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Marg.adhesion", x="",y="")
h_Epith.c.size<-cancerData %>% 
  group_by(Class,Epith.c.size)%>%
  ggplot(aes(x=Epith.c.size, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Epith.c.size", x="",y="")
h_Bare.nuclei<-cancerData %>% 
  group_by(Class,Bare.nuclei)%>%
  ggplot(aes(x=Bare.nuclei, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Bare.nuclei", x="",y="")
h_Bl.cromatin<-cancerData %>% 
  group_by(Class,Bl.cromatin)%>%
  ggplot(aes(x=Bl.cromatin, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Bl.cromatin", x="",y="")
h_Normal.nucleoli<-cancerData %>% 
  group_by(Class,Normal.nucleoli)%>%
  ggplot(aes(x=Normal.nucleoli, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Normal.nucleoli", x="",y="")
h_Mitoses<-cancerData %>% 
  group_by(Class,Mitoses)%>%
  ggplot(aes(x=Mitoses, fill=Class)) +
  geom_histogram(binwidth=0.5)+
  theme(legend.position="bottom")+
  labs(title="Mitoses", x="",y="")

grid.arrange(h_Cl.thickness, h_Cell.size,
             h_Cell.shape,h_Marg.adhesion,
             h_Epith.c.size, h_Bare.nuclei,
             h_Bl.cromatin,h_Normal.nucleoli,
             h_Mitoses,ncol=3)


#4 Box plot
bp_Cl.thickness<-cancerData %>% 
  group_by(Class,Cl.thickness)%>%
  ggplot(aes(x=Class, y=Cl.thickness)) +
  geom_boxplot(outlier.colour="blue", 
               outlier.shape=16,
               outlier.size=2, 
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Cl.thickness", x="",y="")
bp_Cell.size<-cancerData %>% 
  group_by(Class,Cell.size)%>%
  ggplot(aes(x=Class, y=Cell.size)) +
  geom_boxplot(outlier.colour="blue", 
               outlier.shape=16,
               outlier.size=2, 
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Cell.size", x="",y="")
bp_Cell.shape<-cancerData %>% 
  group_by(Class,Cell.shape)%>%
  ggplot(aes(x=Class, y=Cell.shape)) +
  geom_boxplot(outlier.colour="blue", 
               outlier.shape=16,
               outlier.size=2, 
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Cell.shape", x="",y="")
bp_Marg.adhesion<-cancerData %>% 
  group_by(Class,Marg.adhesion)%>%
  ggplot(aes(x=Class, y=Marg.adhesion)) +
  geom_boxplot(outlier.colour="blue",
               outlier.shape=16,
               outlier.size=2, 
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Marg.adhesion", x="",y="")
bp_Epith.c.size<-cancerData %>% 
  group_by(Class,Epith.c.size)%>%
  ggplot(aes(x=Class, y=Epith.c.size)) +
  geom_boxplot(outlier.colour="blue",
               outlier.shape=16,
               outlier.size=2,
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Epith.c.size", x="",y="")
bp_Bare.nuclei<-cancerData %>% 
  group_by(Class,Bare.nuclei)%>%
  ggplot(aes(x=Class, y=Bare.nuclei)) +
  geom_boxplot(outlier.colour="blue",
               outlier.shape=16,
               outlier.size=2,
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Bare.nuclei", x="",y="")
bp_Bl.cromatin<-cancerData %>% 
  group_by(Class,Bl.cromatin)%>%
  ggplot(aes(x=Class, y=Bl.cromatin)) +
  geom_boxplot(outlier.colour="blue",
               outlier.shape=16,
               outlier.size=2,
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Bl.cromatin", x="",y="")
bp_Normal.nucleoli<-cancerData %>% 
  group_by(Class,Normal.nucleoli)%>%
  ggplot(aes(x=Class, y=Normal.nucleoli)) +
  geom_boxplot(outlier.colour="blue", 
               outlier.shape=16,
               outlier.size=2, 
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Normal.nucleoli", x="",y="")
bp_Mitoses<-cancerData %>% 
  group_by(Class,Mitoses)%>%
  ggplot(aes(x=Class, y=Mitoses)) +
  geom_boxplot(outlier.colour="blue",
               outlier.shape=16,
               outlier.size=2,
               notch=FALSE) + 
  theme_minimal() +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=90, 
                                 hjust=0.5, 
                                 vjust=0))+
  labs(title="Mitoses", x="",y="")
grid.arrange(bp_Cl.thickness, bp_Cell.size,
             bp_Cell.shape, bp_Marg.adhesion,
             bp_Epith.c.size, bp_Bare.nuclei,
             bp_Bl.cromatin, bp_Normal.nucleoli,
             bp_Mitoses, ncol=3)



#Normalisasi
normalization <- function(x){
  (x - min(x)) / (max(x)-min(x))
}

cancerData <- cancerData %>%
  mutate(Cl.thickness = normalization(Cl.thickness), 
         Cell.size = normalization(Cell.size), 
         Cell.shape = normalization(Cell.shape), 
         Marg.adhesion = normalization(Marg.adhesion), 
         
         Epith.c.size = normalization(Epith.c.size), 
         Bare.nuclei = normalization(Bare.nuclei), 
         Bl.cromatin = normalization(Bl.cromatin),
         
         Normal.nucleoli = normalization(Normal.nucleoli), 
         Mitoses = normalization(Mitoses), 
         Class = as.numeric(Class)-1)
head(cancerData)

#Type
typeof(cancerData$Cl.thickness)
typeof(cancerData$Cell.size)
typeof(cancerData$Cell.shape)
typeof(cancerData$Marg.adhesion)
typeof(cancerData$Epith.c.size)
typeof(cancerData$Bare.nuclei)
typeof(cancerData$Bl.cromatin)
typeof(cancerData$Normal.nucleoli)
typeof(cancerData$Mitoses)
typeof(cancerData$Class)

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

#Model Matrix
str(cancerTrainingData)
model <- model.matrix( ~ Class + Cl.thickness + Cell.size + Cell.shape + Marg.adhesion +
    Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
  data = cancerTrainingData )
summary(model)
head(model)

#Model NN
nn_model <- neuralnet(Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion +
                             Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
                           data=model, hidden=c(2,1), linear.output=FALSE, threshold=0.01)
nn_model$result.matrix
plot(nn_model)


#Predict: Test the resulting output
test_result <- subset(cancerTestData, select = c("Cl.thickness", "Cell.size", "Cell.shape", "Marg.adhesion",
                                                 "Epith.c.size", "Bare.nuclei", "Bl.cromatin", "Normal.nucleoli", "Mitoses"))
head(test_result)
nn.results <- compute(nn_model, test_result)
results <- data.frame(actual = cancerTestData$Class, prediction = nn.results$net.result)

#Predict: Confusion Matrix
roundedresults<-sapply(results,round,digits=0)
roundedresults_dataframe=data.frame(roundedresults)
attach(roundedresults_dataframe)
table(actual,prediction)

#Measure Performance
confusion_matrix=table(actual,prediction)
TP <- confusion_matrix[1,1]
FP <- confusion_matrix[1,2]
FN <- confusion_matrix[2,1]
TN <- confusion_matrix[2,2]

#1. Accuracy
accuracy <- function(TP,FP,FN,TN){
  accuracy=(TP+TN)/(TP+FP+FN+TN)
  return(accuracy)
}
accuracy(TP,FP,FN,TN)

#2. Classification Error Rate
CER <- function(TP,FP,FN,TN){
  CER=(FP+FN)/(TP+FP+FN+TN)
  return(CER)
}
CER(TP,FP,FN,TN)

#3. Precision
precision <- function(TP,FP,FN,TN){
  precision=(TP)/(TP+FP)
  return(precision)
}
precision(TP,FP,FN,TN)

#4. Sensitivity
sensitivity <- function(TP,FP,FN,TN){
  sensitivity=(TP)/(TP+FN)
  return(sensitivity)
}
sensitivity(TP,FP,FN,TN)

#5. Specificity
specificity <- function(TP,FP,FN,TN){
  specificity=(TN)/(TN+FP)
  return(specificity)
}
specificity(TP,FP,FN,TN)

#6. F1
F1 <- function(precision,sensitivity){
  F1=(2*precision*sensitivity)/(precision+sensitivity)
  return(F1)
}
F1(1,0.96)




cancerDataMutate <- cancerData %>%
  mutate(Class = as.integer(Class) - 1, 
         Class = ifelse(Class == 1, TRUE, FALSE))



#Tune Hyper-parameters
set.seed(489)
# 2-Hidden Layers, Layer-1 2-neurons, Layer-2, 1-neuron
nn_model2 <- neuralnet(Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion +
                      Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
                     data = cancerDataMutate, 
                     linear.output = FALSE, 
                     err.fct = 'ce', 
                     likelihood = 
                       TRUE, hidden = c(2,1))

# 2-Hidden Layers, Layer-1 2-neurons, Layer-2, 2-neurons
set.seed(489)
nn_model3 <- nn_model2 <- neuralnet(Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion +
                                  Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
                                data = cancerDataMutate, 
                                linear.output = FALSE, 
                                err.fct = 'ce', 
                                likelihood = TRUE, 
                                hidden = c(2,2))

# 2-Hidden Layers, Layer-1 1-neuron, Layer-2, 2-neuron
set.seed(489)
nn_model4 <- nn_model2 <- neuralnet(Class ~ Cl.thickness + Cell.size + Cell.shape + Marg.adhesion +
                                  Epith.c.size + Bare.nuclei + Bl.cromatin + Normal.nucleoli + Mitoses, 
                                data = cancerDataMutate, 
                                linear.output = FALSE, 
                                err.fct = 'ce', 
                                likelihood = TRUE, 
                                hidden = c(1,2))

# Bar plot of results
class_nn <- tibble('Network' = rep(c("NN1", "NN2", "NN3", "NN4"), each = 3), 
                       'Metric' = rep(c('AIC', 'BIC', 'ce Error * 100'), length.out = 12),
                       'Value' = c(nn_model$result.matrix[4,1], nn_model$result.matrix[5,1], 
                                   100*nn_model$result.matrix[1,1], nn_model2$result.matrix[4,1], 
                                   nn_model2$result.matrix[5,1], 100*nn_model2$result.matrix[1,1],
                                   nn_model3$result.matrix[4,1], nn_model3$result.matrix[5,1], 
                                   100*nn_model3$result.matrix[1,1], nn_model4$result.matrix[4,1], 
                                   nn_model4$result.matrix[5,1], 100*nn_model4$result.matrix[1,1]))

class_nn %>%
  ggplot(aes(Network, Value, fill = Metric)) +
  geom_col(position = 'dodge')  +
  ggtitle("AIC, BIC, and Cross-Entropy Error of the Classification ANNs", "Note: ce Error displayed is 100 times its true value")
