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
  d2$gpd_dryweight_cor=d2$growth_per_day-d2$InitialSize
  d2$gpd_dryweight_cor # A super simple way of correcting has earlier shown great predictability
  
  fit <- lmer(growth_per_day ~ InitialSize + (1|Clover), data=d2) 
  ycorr <- d2$growth_per_day - model.matrix( ~ InitialSize, data=d2) %*% fixef(fit)
  d2$gpd_dryweight_cor <- ycorr #this is the new corrected dry weight
  }

# Take out round 1, replicate 2
d2$roundRep <- paste(d2$Round, d2$Replicate, sep='_')


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

# Make the matrices that goes into the model
{
  CloverDesign <- model.matrix(~0+d6$Clovershort)
  GcloverReps <- CloverDesign %*% GRM1 %*% t(CloverDesign) 
  #GcloverReps is a G-matrix that will match the size and lay-out in the data and can be used in BGLR at the K=
  
  
  #Add a clover effect to capture non-additve variance that makes broad sense heritability.
  #CloverIndep catch clover-effects without relationships (clovers are independent), that is another matrix that can go into the model. 
  CloverIndep <- CloverDesign %*% t(CloverDesign)
  dim(CloverIndep)
}

Correctedforallfixed <- lmer(gpd_dryweight_cor ~ factor(roundRep) + factor(NS) + factor(EW) + factor(Rhizobium) + (1|Clover), data=d6) 
summary(Correctedforallfixed)

ycorr <- d6$gpd_dryweight_cor - model.matrix( ~ factor(roundRep) + factor(NS) + factor(EW)  + factor(Rhizobium), data=d6) %*% fixef(Correctedforallfixed)
d6$CorrectedPheno <- ycorr  

GP_GBLUP<-function(testpop){
  yNA=y
  yNA[testpop]=NA
  fixedmod=model.matrix(~factor(d6$roundRep)+factor(d6$NS)+factor(d6$EW)+factor(d6$Rhizobium))
  ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
  GBLUP=BGLR(y=yNA,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("GBLUP",round))
  matrix=cbind(as.character(d6$Clovershort[testpop]),as.numeric(d6$CorrectedPheno[testpop]),as.numeric(GBLUP$ETA[[1]]$u[testpop]),as.numeric(GBLUP$ETA[[2]]$u[testpop]),(as.numeric(GBLUP$ETA[[1]]$u[testpop])+as.numeric(GBLUP$ETA[[2]]$u[testpop])))
  colnames(matrix)=c("ID", "Observed", "1st GEBV contribution","2nd GEBV contribution (clover independent)","totalGEBV")
  return(matrix)
}


#Load groups from file

f=read.table(args[3],fill = TRUE)

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

testpop1_idx=which(d6$Clovershort %in% testpop1)
testpop2_idx=which(d6$Clovershort %in% testpop2)
testpop3_idx=which(d6$Clovershort %in% testpop3)
testpop4_idx=which(d6$Clovershort %in% testpop4)
testpop5_idx=which(d6$Clovershort %in% testpop5)
testpop6_idx=which(d6$Clovershort %in% testpop6)

tests=list(testpop1_idx,testpop2_idx,testpop3_idx,testpop4_idx,testpop5_idx,testpop6_idx)


#  Calculate heritability
print("Calculating heritability")
y=d6$gpd_dryweight_cor
fixedmod=model.matrix(~factor(d6$roundRep)+factor(d6$NS)+factor(d6$EW)+factor(d6$Rhizobium))
ETA=list(list(K=GcloverReps,model="RKHS"),list(K=CloverIndep,model="RKHS"),list(X=fixedmod,model="FIXED")) 
GBLUP=BGLR(y=y,response_type = "gaussian",ETA=ETA,nIter=20000,burnIn = 5000,verbose=F,saveAt=paste("heritabilityestimation",round))

h2=(GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)/(GBLUP$varE+GBLUP$ETA[[1]]$varU+GBLUP$ETA[[2]]$varU)
write.table(h2,paste("h2",round,".txt",sep=""),sep="\t")


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

filename=paste("Correlation_GBLUP_GPD",round,".txt",sep="")
write.table(correlation,filename,sep="\t")

filename1=paste("Predictions_GBLUP_GPD",round,".txt",sep="")
write.table(Combined,filename1,sep="\t")


  