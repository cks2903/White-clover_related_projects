####################################################################################
# Genomic prediction fitting fixed effects on yield data with replicates.          #
####################################################################################
#                       GBLUP                                                      #
####################################################################################

#setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GP_based_on_RNASNPs/BacktoBasicYield_092019_GP/")
# this is a script to run BayesB on replicate data with predefined division into training and testing groups
library(lme4)
library(BGLR)
library("parallel")
library("methods")
library("Matrix")
args=commandArgs(trailingOnly = TRUE)
print(args)
round=args[2]

# read in data
{
  d <- read.csv("/home/cks/NChain/faststorage/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/greenhouse_area.csv", header = TRUE, sep = ",")
  head(d)
  
  f=read.csv("/home/cks/NChain/faststorage/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/2018_weight.csv",header=T,sep=";")
  colnames(f)[1]="Barcode"
  head(f)
  
  df=merge(d,f,by="Barcode")
  head(df)
  d=df
  
}

#Calculate growth pr. day
{
  d$days_of_growth <- as.Date(d$harvest_date, format = "%d/%m/%y") - as.Date(d$inoculation_date, format = "%d/%m/%y")
  d$growth_per_day <- d$weight.y/as.numeric(d$days_of_growth)
}


# Sort out outliers for initialsize, plants that don't grow, plants that were inoculated with no rhizobium, SM73 or had n_stolon=0. These are errors and don't grow
{
  
  d0=na.omit(d)
  hist(d0$InitialSize)
  d1 <- subset(d0, InitialSize<40000)
  hist(d1$InitialSize)
  
  hist(d1$growth_per_day,breaks = 500)
  
  sortout=round(length(d1$growth_per_day)*0.05,0)
  sortout #remove 175 closest to 0 in growth pr. day
  idxout=order(d1$growth_per_day,decreasing=F)[1:sortout]
  d2=d1[-idxout,]
  hist(d2$growth_per_day,breaks = 500)
  
  d2=d2[-which(d2$rhizobium=="SM73"),]
  d2=d2[-which(d2$rhizobium=="NO"),]
  d2=d2[-which(d2$n_stolons==0),]
  d2=d2[-which(d2$n_stolons==1),]
  d2=d2[-which(d2$n_stolons==2),]
  d2=d2[-which(d2$n_stolons==3),]
}

# Correct GPD for the part of initial size that does not carry genetic information
{
  lm.fit <- lm(d2$InitialSize ~ d2$Clover)
  summary(lm.fit)
  d2$residuals <- round(lm.fit$residuals,2)
  
  d2$GPD_initialsize_cor=d2$growth_per_day-d2$residuals
  head(d2)
  d2$GPD_initialsize_cor
  
  fit <- lmer(growth_per_day ~ residuals + (1|Clover), data=d2) 
  ycorr <- d2$growth_per_day - model.matrix( ~ residuals, data=d2) %*% fixef(fit)
  d2$gpd_dryweight_cor <- ycorr #this is the new corrected dry weight
}

# Make roundrep variable
{
d2$roundRep <- paste(d2$Round, d2$Replicate, sep='_')
}

# Fix data matrix so that names are the same as in GRM, that is 118 individuals in each
{
  GRM=read.table(args[1],sep=",",header=T)
  dim(GRM)
  d2$Clovershort <- strtrim(d2$Clover,8)
  d3=d2[order(d2$Clovershort,decreasing=F),]
  length(unique(d3$Clovershort)) #149
  d4=d3[which(d3$Clovershort %in% colnames(GRM)),]
  length(unique(d4$Clovershort)) #112 unique genotypes with GPD data
  remove=GRM[which(colnames(GRM) %in% d4$Clovershort),which(colnames(GRM) %in% d4$Clovershort)]
  print(remove)
  GRM1=GRM[which(colnames(GRM) %in% d4$Clovershort),which(colnames(GRM) %in% d4$Clovershort)]
  dim(GRM1)
  nrow(GRM1)==length(unique(d4$Clovershort))
  GRM1=data.matrix(GRM1)
}

#Aberpearl_07 contribute 700 of the datapoints and thus influence the variance a lot. Cut down Aberpearl_07 data so that we have only 6 different Rhizobia left like the other clovers
{
Aberpearl_07=which(d4$Clovershort=="Aearl_07")
Inocolums=unique(d4$Rhizobium[Aberpearl_07])
set.seed(15)
sample=sample(Inocolums,6)

which(d4$Rhizobium[Aberpearl_07]==sample[1]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[2]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[3]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[4]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[5]) #4
which(d4$Rhizobium[Aberpearl_07]==sample[6]) #4

remove=which((d4$Rhizobium[Aberpearl_07] %in% sample)==FALSE)

d4=d4[-Aberpearl_07[remove],]
nrow(d4)


#Remove round 1 replicate 2
d6=d4[-which(d4$roundRep=="1_2"),]
nrow(d6)


# Sort out genotypes that do not have at least 10 replicates as they will disturb the analysis
Genotypes=(unique(d6$Clovershort))

for (genotype in Genotypes){
  idx=which(d6$Clovershort==genotype)
  if (length(idx)<10){
    d6=d6[-idx,]
  }
}


# Remove genotypes not included from GRM
{
  GRM1=GRM1[which(colnames(GRM1) %in% d6$Clovershort),which(colnames(GRM1) %in% d6$Clovershort)]
  nrow(GRM1)==length(unique(d6$Clovershort))
  GRM1=data.matrix(GRM1)
}


# Divide into 6 populations for GP
set.seed(NULL)
tst=sample(1:length(unique(d6$Clovershort)),size=length(unique(d6$Clovershort)),replace=FALSE) 
k=6
testing_pop=split(tst, sort(tst%%k))
  
tst1=testing_pop[1]$'0'
tst2=testing_pop[2]$'1'
tst3=testing_pop[3]$'2'
tst4=testing_pop[4]$'3'
tst5=testing_pop[5]$'4'
tst6=testing_pop[6]$'5'
  
testpop1=unique(d6$Clovershort)[tst1]
testpop2=unique(d6$Clovershort)[tst2]
testpop3=unique(d6$Clovershort)[tst3]
testpop4=unique(d6$Clovershort)[tst4]
testpop5=unique(d6$Clovershort)[tst5]
testpop6=unique(d6$Clovershort)[tst6]
  
grouping=list(testpop1,testpop2,testpop3,testpop4,testpop5,testpop6)
name=paste("grouping",round,".txt",sep="")
sink(name)
print(grouping)
sink()


############################################################



# Now remove replicates so each genotype has a maximum of of desired number (maxreplicates) 

removereplicates <- function(maxreplicates,dataframe){
   for (genotype in Genotypes){
    replicateidx=which(dataframe$Clovershort==genotype)
    if (length(replicateidx)>maxreplicates){
      numbertoremove=length(replicateidx)-maxreplicates
      remove=sample(replicateidx,numbertoremove)
      dataframe=dataframe[-remove,]
    }
  }
  return(dataframe)
}


# Function that takes in a dataframe, prepare genomic prediction matrices, and do GP
GP <- function(maxreplicates,dataframe){

  # Make the matrices that goes into the model
  CloverDesign <- model.matrix(~0+dataframe$Clovershort)
  GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
  #GcloverReps is a G-matrix that will match the size and lay-out in the data and can be used in BGLR at the K=

  #Add a clover effect to capture non-additve variance that makes broad sense heritability.
  #CloverIndep catch clover-effects without relationships (clovers are independent), that is another matrix that can go into the model. 
  CloverIndep <- CloverDesign %*% t(CloverDesign)
  dim(CloverIndep)
  dataframe$Rhizobium <- factor(dataframe$Rhizobium)
  dataframe$Clover <- factor(dataframe$Clover)
  dataframe$NS <- factor(dataframe$NS)
  dataframe$EW <- factor(dataframe$EW)

  Correctedforallfixed <- lmer(gpd_dryweight_cor ~ factor(roundRep) + factor(NS) + factor(EW) + factor(Rhizobium) + (1|Clover), data=dataframe) 
  summary(Correctedforallfixed)

  ycorr <- dataframe$gpd_dryweight_cor - model.matrix( ~ factor(roundRep) + factor(NS) + factor(EW)  + factor(Rhizobium), data=dataframe) %*% fixef(Correctedforallfixed)
  dataframe$CorrectedPheno <- ycorr  

  GP_GBLUP<-function(testpop){
    yNA=y
    yNA[testpop]=NA
    fixedmod=model.matrix(~factor(dataframe$roundRep)+factor(dataframe$NS)+factor(dataframe$EW)+factor(dataframe$Rhizobium))
    ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
    GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
    matrix=cbind(as.character(dataframe$Clovershort[testpop]),as.numeric(dataframe$CorrectedPheno[testpop]),as.numeric(GBLUP$ETA[[1]]$u[testpop]),as.numeric(GBLUP$ETA[[2]]$u[testpop]),(as.numeric(GBLUP$ETA[[1]]$u[testpop])+as.numeric(GBLUP$ETA[[2]]$u[testpop])))
    colnames(matrix)=c("ID", "Observed", "1st GEBV contribution","2nd GEBV contribution (clover independent)","totalGEBV")
    return(matrix)
  }



#Find indexes for test population  
  testpop1_idx=which(dataframe$Clovershort %in% testpop1)
  testpop2_idx=which(dataframe$Clovershort %in% testpop2)
  testpop3_idx=which(dataframe$Clovershort %in% testpop3)
  testpop4_idx=which(dataframe$Clovershort %in% testpop4)
  testpop5_idx=which(dataframe$Clovershort %in% testpop5)
  testpop6_idx=which(dataframe$Clovershort %in% testpop6)

  
  tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
  

#  Calculate heritability
  print("Calculating heritability")
  y=dataframe$gpd_dryweight_cor
  fixedmod=model.matrix(~factor(dataframe$roundRep)+factor(dataframe$NS)+factor(dataframe$EW)+factor(dataframe$Rhizobium))
  ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation",round))

  h2=(GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
  write.table(h2,paste("h2",round,maxreplicates,".txt",sep=""),sep="\t")


  # Prediction
  print("Starting GBLUP prediction")
  results=mclapply(tests,GP_GBLUP)

  first=results[[1]]
  second=results[[2]]
  third=results[[3]]
  fourth=results[[4]]
  fifth=results[[5]]
  sixth=results[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)

  All1=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)
  All4=aggregate(as.numeric(All[, 5]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All3=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  correlation=cor(All1[,2],All4[,2]) #means of replicates
  correlation

  Combined=cbind(All1[,1],All1[,2],All2[,2],All3[,2],All4[,2])
  colnames(Combined)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs Cloverreps","GEBVs CloverIndep","Total GEBVs")

  filename=paste("Correlation_GBLUP_GPD",round,maxreplicates,".txt",sep="")
  write.table(correlation,filename,sep="\t")

  filename1=paste("Predictions_GBLUP_GPD",round,maxreplicates,".txt",sep="")
  write.table(Combined,filename1,sep="\t")

return()
}


GP_GBLUP<-function(testpop){
  yNA=y
  yNA[testpop]=NA
  fixedmod=model.matrix(~factor(dataframe$roundRep)+factor(dataframe$NS)+factor(dataframe$Rhizobium))
  ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(dataframe$Clovershort[testpop]),as.numeric(dataframe$CorrectedPheno[testpop]),as.numeric(GBLUP$ETA[[1]]$u[testpop]),as.numeric(GBLUP$ETA[[2]]$u[testpop]),(as.numeric(GBLUP$ETA[[1]]$u[testpop])+as.numeric(GBLUP$ETA[[2]]$u[testpop])))    colnames(matrix)=c("ID", "Observed", "1st GEBV contribution","2nd GEBV contribution (clover independent)","totalGEBV")
  return(matrix)
  }



#Find indexes for test population  
testpop1_idx=which(dataframe$Clovershort %in% testpop1)
testpop2_idx=which(dataframe$Clovershort %in% testpop2)
testpop3_idx=which(dataframe$Clovershort %in% testpop3)
testpop4_idx=which(dataframe$Clovershort %in% testpop4)
testpop5_idx=which(dataframe$Clovershort %in% testpop5)
testpop6_idx=which(dataframe$Clovershort %in% testpop6)

  
tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)
  

#  Calculate heritability
  print("Calculating heritability")
  y=dataframe$gpd_dryweight_cor
  fixedmod=model.matrix(~factor(dataframe$roundRep)+factor(dataframe$NS)+factor(dataframe$Rhizobium))
  ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation",round))

  h2=(GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
  write.table(h2,paste("h2",round,maxreplicates,".txt",sep=""),sep="\t")


  # Prediction
  print("Starting GBLUP prediction")
  results=mclapply(tests,GP_GBLUP)

  first=results[[1]]
  second=results[[2]]
  third=results[[3]]
  fourth=results[[4]]
  fifth=results[[5]]
  sixth=results[[6]]
  All=rbind(first,second,third,fourth,fifth,sixth)

  All1=aggregate(as.numeric(All[, 2]), list(All[,1]), mean)
  All4=aggregate(as.numeric(All[, 5]), list(All[,1]), mean)
  All2=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  All3=aggregate(as.numeric(All[, 3]), list(All[,1]), mean)
  correlation=cor(All1[,2],All4[,2]) #means of replicates
  correlation

  Combined=cbind(All1[,1],All1[,2],All2[,2],All3[,2],All4[,2])
  colnames(Combined)=c("Individual","ObservedCorrectedforFixedeffects","GEBVs Cloverreps","GEBVs CloverIndep","Total GEBVs")

  filename=paste("Correlation_GBLUP_GPD",round,maxreplicates,".txt",sep="")
  write.table(correlation,filename,sep="\t")

  filename1=paste("Predictions_GBLUP_GPD",round,maxreplicates,".txt",sep="")
  write.table(Combined,filename1,sep="\t")

return()
}






#Apply so maximum of replicates is 10
{
Only10reps=removereplicates(10,d6)
Only10reps=as.data.frame(Only10reps)
nrow(Only10reps)
GP(10,Only10reps)
}

#Apply so maximum is 9
{
Only9reps=removereplicates(9,Only10reps)
Only9reps=as.data.frame(Only9reps)
nrow(Only9reps)
GP(9,Only9reps)
}

#Apply so maximum of replicates is 8
{
Only8reps=removereplicates(8,Only9reps)
Only8reps=as.data.frame(Only8reps)
nrow(Only8reps)
GP(8,Only8reps)
}

#Apply so maximum of replicates is 7
{
Only7reps=removereplicates(7,Only8reps)
Only7reps=as.data.frame(Only7reps)
nrow(Only7reps)
GP(7,Only7reps)
}

#Apply so maximum of replicates is 6
{
Only6reps=removereplicates(6,Only7reps)
Only6reps=as.data.frame(Only6reps)
nrow(Only6reps)
GP(6,Only6reps)
}

#Apply so maximum of replicates is 5
{
Only5reps=removereplicates(5,Only6reps)
Only5reps=as.data.frame(Only5reps)
nrow(Only5reps)
GP(5,Only5reps)
}

#Apply so maximum of replicates is 4
{
Only4reps=removereplicates(4,Only5reps)
Only4reps=as.data.frame(Only4reps)
nrow(Only4reps)
GP(4,Only4reps)
}

#Apply so maximum of replicates is 3
{
Only3reps=removereplicates(3,Only4reps)
Only3reps=as.data.frame(Only3reps)
nrow(Only3reps)
GP(3,Only3reps)
}

#Apply so maximum of replicates is 2 is not possible as "fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients"

