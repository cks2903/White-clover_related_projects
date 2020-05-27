#########################################################################
#                 Using machine learning for genomic prediction         #
#                                                                       #
#########################################################################

#install.packages("caret")
#install.packages("e1071",dependencies = T)
library(caret)
args=commandArgs(trailingOnly = TRUE)

##An example comes with the mlbench package
{
library(mlbench)
#data("iris") #5 columns sepal length, sepal width. petal length, petal widt and species (classification)
#iris$Species = as.factor(iris$Species)
#To divide data into training and testing (25% of the data)
#inTraining <- createDataPartition(iris$Species, p = .75, list = FALSE)
#training <- iris[ inTraining,]
#testing  <- iris[-inTraining,]


#Create the Cross validation method 
#FitControl=trainControl(method="repeatedcv",number=10,repeats=10)		# A 10-fold CV

# Train a model
#set.seed(12345)


#model=train(training[,2:ncol(training)],as.factor(training[,1]),method= "ranger",importance=permutation) #fit training
#model

#gbmFit1 <- train(Species ~ ., data = training, 
                 #method = "gbm", 
                 #trControl = FitControl,
                 #verbose = FALSE)

#gbmFit1

#Accuracy is used by Caret to chose the final model. It is the overall agreement rate btw. the CV methods
#Kappat is another statistical method used for accessing models with categorical variables 
#CARET chose the first model with an interaction depth of 1, number of trees at 50, an accuracy of 97% and a Kappa of 95%.
#The models very in shrinkage 


#Finally we can use the training model to predict classifications and probabilities for the test data set
#predictions<-predict(object=model,testing[,2:ncol(testing)])
#predictions
}

#######################try random forest model on your data######################################

#Load phenotypes, corrected for fixed effects and GRM (outliers filtered out and LD filtered)
{
#setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GP_based_on_RNASNPs/BacktoBasicYield_092019_GP")  

GRM=read.table("GRM_outliersremoved_RNAseq_LDfiltered_20190902.csv",sep=";",header=T)  
dim(GRM)

pheno=read.table("CorrectedGPDPhenotypes_20191009.csv",sep=",",header=T)
dim(pheno)

#remove individuals in GRM not in pheno
pheno[,1]

NAMES<- strtrim(pheno$Accession,8)
pheno$Accession=NAMES

Remove=which(colnames(GRM) %in% NAMES==FALSE)
GRM1=GRM[-Remove,-Remove]
dim(GRM1)

GRM1=data.matrix(GRM1)
}

#then make GRM 1292x1292 in dimension
{
  CloverDesign <- model.matrix(~0+pheno$Accession)
  GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
  #GcloverReps is a G-matrix that will match the size and lay-out in the data and can be used in BGLR at the K=
}

#Combine genotypes and phenotypes into one large dataframe
{
cloverdata=cbind(pheno,GcloverReps)
dim(cloverdata)    
}

# Divide into predefined populations used for 6-fold CV
{
f=read.table(args[1],fill = TRUE)

linewherenewgroupcomes=vector()
for (i in 1:nrow(f)){
  beginning=as.character(f[i,1])
  secondchr=strsplit(beginning,"")[[1]][2]
  if (secondchr=="["){
    linewherenewgroupcomes=append(linewherenewgroupcomes,i)
  }
}

group1=f[(linewherenewgroupcomes[1]+1):(linewherenewgroupcomes[1+1]-1),2:ncol(f)]
testpop1=c(unique(as.character(group1[,1])),unique(as.character(group1[,2])),unique(as.character(group1[,3])),unique(as.character(group1[,4])),unique(as.character(group1[,5])),unique(as.character(group1[,6])))

group2=f[(linewherenewgroupcomes[2]+1):(linewherenewgroupcomes[2+1]-1),2:ncol(f)]
testpop2=c(unique(as.character(group2[,1])),unique(as.character(group2[,2])),unique(as.character(group2[,3])),unique(as.character(group2[,4])),unique(as.character(group2[,5])),unique(as.character(group2[,6])))

group3=f[(linewherenewgroupcomes[3]+1):(linewherenewgroupcomes[3+1]-1),2:ncol(f)]
testpop3=c(unique(as.character(group3[,1])),unique(as.character(group3[,2])),unique(as.character(group3[,3])),unique(as.character(group3[,4])),unique(as.character(group3[,5])),unique(as.character(group3[,6])))

group4=f[(linewherenewgroupcomes[4]+1):(linewherenewgroupcomes[4+1]-1),2:ncol(f)]
testpop4=c(unique(as.character(group4[,1])),unique(as.character(group4[,2])),unique(as.character(group4[,3])),unique(as.character(group4[,4])),unique(as.character(group4[,5])),unique(as.character(group4[,6])))

group5=f[(linewherenewgroupcomes[5]+1):(linewherenewgroupcomes[5+1]-1),2:ncol(f)]
testpop5=c(unique(as.character(group5[,1])),unique(as.character(group5[,2])),unique(as.character(group5[,3])),unique(as.character(group5[,4])),unique(as.character(group5[,5])),unique(as.character(group5[,6])))

group6=f[(linewherenewgroupcomes[6]+1):nrow(f),2:ncol(f)]
testpop6=c(unique(as.character(group6[,1])),unique(as.character(group6[,2])),unique(as.character(group6[,3])),unique(as.character(group6[,4])),unique(as.character(group6[,5])),unique(as.character(group6[,6])))

testpop1_idx=which(pheno$Accession %in% testpop1)
testpop2_idx=which(pheno$Accession %in% testpop2)
testpop3_idx=which(pheno$Accession %in% testpop3)
testpop4_idx=which(pheno$Accession %in% testpop4)
testpop5_idx=which(pheno$Accession %in% testpop5)
testpop6_idx=which(pheno$Accession %in% testpop6)

tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
}


# Machine learning using the first group as testing and group2,3,4,5,6 as training
#cloverdata$GPD_corrected = as.factor(cloverdata$GPD_corrected)

training <- cloverdata[-testpop1_idx,]
dim(training)
testing  <- cloverdata[testpop1_idx,]
dim(testing)


model=train(training[,3:ncol(training)],training[,2],method= "ranger") #fit training
model

predictions_group1 =predict(object=model, testing[,3:ncol(testing)])
predictions_group1

Predicted_obs1=cbind(pheno[testpop1_idx,1],predictions_group1,pheno[testpop1_idx,2])
colnames(Predicted_obs1)=c("Individual","Prediction","Corrected_pheno")


MeanPred=aggregate(as.numeric(Predicted_obs1[, 2]), list(Predicted_obs1[,1]), mean)
MeanObs=aggregate(as.numeric(Predicted_obs1[, 3]), list(Predicted_obs1[,1]), mean)

Predicted_obs1_1=cbind(MeanPred,MeanObs[,2])
colnames(Predicted_obs1_1)=c("Individual","Prediction","Corrected_pheno")


# Machine learning using the second group as testing and group1,3,4,5,6 as training
training <- cloverdata[-testpop2_idx,]
dim(training)
testing  <- cloverdata[testpop2_idx,]
dim(testing)



model=train(training[,3:ncol(training)],training[,2],method= "ranger") #fit training
model

predictions_group2 =predict(object=model, testing[,3:ncol(testing)])
predictions_group2


Predicted_obs2=cbind(pheno[testpop2_idx,1],predictions_group2,pheno[testpop2_idx,2])
colnames(Predicted_obs2)=c("Individual","Prediction","Corrected_pheno")


MeanPred=aggregate(as.numeric(Predicted_obs2[, 2]), list(Predicted_obs2[,1]), mean)
MeanObs=aggregate(as.numeric(Predicted_obs2[, 3]), list(Predicted_obs2[,1]), mean)

Predicted_obs2_2=cbind(MeanPred,MeanObs[,2])
colnames(Predicted_obs2_2)=c("Individual","Prediction","Corrected_pheno")


# Machine learning using the third group as testing and group1,2,4,5,6 as training
training <- cloverdata[-testpop3_idx,]
dim(training)
testing  <- cloverdata[testpop3_idx,]
dim(testing)


model=train(training[,3:ncol(training)],training[,2],method= "ranger") #fit training
model

predictions_group3 =predict(object=model, testing[,3:ncol(testing)])
predictions_group3

Predicted_obs3=cbind(pheno[testpop3_idx,1],predictions_group3,pheno[testpop3_idx,2])
colnames(Predicted_obs3)=c("Individual","Prediction","Corrected_pheno")


MeanPred=aggregate(as.numeric(Predicted_obs3[, 2]), list(Predicted_obs3[,1]), mean)
MeanObs=aggregate(as.numeric(Predicted_obs3[, 3]), list(Predicted_obs3[,1]), mean)

Predicted_obs3_3=cbind(MeanPred,MeanObs[,2])
colnames(Predicted_obs3_3)=c("Individual","Prediction","Corrected_pheno")





# Machine learning using the fourth group as testing and group1,2,3,5,6 as training
training <- cloverdata[-testpop4_idx,]
dim(training)
testing  <- cloverdata[testpop4_idx,]
dim(testing)


model=train(training[,3:ncol(training)],training[,2],method= "ranger") #fit training
model

predictions_group4 =predict(object=model, testing[,3:ncol(testing)])
predictions_group4

Predicted_obs4=cbind(pheno[testpop4_idx,1],predictions_group4,pheno[testpop4_idx,2])
colnames(Predicted_obs4)=c("Individual","Prediction","Corrected_pheno")


MeanPred=aggregate(as.numeric(Predicted_obs4[, 2]), list(Predicted_obs4[,1]), mean)
MeanObs=aggregate(as.numeric(Predicted_obs4[, 3]), list(Predicted_obs4[,1]), mean)

Predicted_obs4_4=cbind(MeanPred,MeanObs[,2])
colnames(Predicted_obs4_4)=c("Individual","Prediction","Corrected_pheno")


# Machine learning using the fifth group as testing and group1,2,3,4,6 as training
training <- cloverdata[-testpop5_idx,]
dim(training)
testing  <- cloverdata[testpop5_idx,]
dim(testing)


model=train(training[,3:ncol(training)],training[,2],method= "ranger") #fit training
model

predictions_group5 =predict(object=model, testing[,3:ncol(testing)])
predictions_group5

Predicted_obs5=cbind(pheno[testpop5_idx,1],predictions_group5,pheno[testpop5_idx,2])
colnames(Predicted_obs5)=c("Individual","Prediction","Corrected_pheno")


MeanPred=aggregate(as.numeric(Predicted_obs5[, 2]), list(Predicted_obs5[,1]), mean)
MeanObs=aggregate(as.numeric(Predicted_obs5[, 3]), list(Predicted_obs5[,1]), mean)

Predicted_obs5_5=cbind(MeanPred,MeanObs[,2])
colnames(Predicted_obs5_5)=c("Individual","Prediction","Corrected_pheno")




# Machine learning using the sixth group as testing and group1,2,3,4,5 as training
training <- cloverdata[-testpop6_idx,]
dim(training)
testing  <- cloverdata[testpop6_idx,]
dim(testing)


model=train(training[,3:ncol(training)],training[,2],method= "ranger") #fit training
model

predictions_group6 =predict(object=model, testing[,3:ncol(testing)])
predictions_group6

Predicted_obs6=cbind(pheno[testpop6_idx,1],predictions_group6,pheno[testpop6_idx,2])
colnames(Predicted_obs6)=c("Individual","Prediction","Corrected_pheno")


MeanPred=aggregate(as.numeric(Predicted_obs6[, 2]), list(Predicted_obs6[,1]), mean)
MeanObs=aggregate(as.numeric(Predicted_obs6[, 3]), list(Predicted_obs6[,1]), mean)

Predicted_obs6_6=cbind(MeanPred,MeanObs[,2])
colnames(Predicted_obs6_6)=c("Individual","Prediction","Corrected_pheno")


#Combine
Combined=rbind(Predicted_obs1_1,Predicted_obs2_2,Predicted_obs3_3,Predicted_obs4_4,Predicted_obs5_5,Predicted_obs6_6)
Correlation=cor(Combined[,2],Combined[,3])

filename=args[2]

write.table(Combined,filename,sep="\t",row.names=F,col.names=T,quote=F)
