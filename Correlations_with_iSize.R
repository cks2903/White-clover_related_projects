#####################################################
#       Correlations between iSize and GPD traits
#####################################################

library(corrgram)
library(lme4)
library(BGLR)
library("parallel")
library("methods")
library("Matrix")
library(agricolae)

setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GWAS results/GWAS_20200206")

# read in data
{
  d <- read.csv("greenhouse_area.csv", header = TRUE, sep = ",")
  head(d)
  
  f=read.csv("2018_weight.csv",header=T,sep=";")
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

# Take out round 1, replicate 2
d2$roundRep <- paste(d2$Round, d2$Replicate, sep='_')


# Fix data matrix so that names are the same as in GRM, that is 118 individuals in each
{
  GRM=read.table("GRM_LDfilteredwithingene_20191211_nooutliers.csv",sep=",",header=T)
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

#Remove round 1 replicate 2
d6=d4[-which(d4$roundRep=="1_2"),]
nrow(d6)


nrow(d6)/length(unique(d6$Clovershort)) #15 observations pr. genotype on average. But Aberpearl_07 has 405 of these
(nrow(d6)-405)/(length(unique(d6$Clovershort))-1)  # 11.5 observations pr. genotype that is not Aberpearl_07


 # Correcting for initial size as fixed effect without only using residual effects
{
  fit <- lmer(growth_per_day ~ InitialSize + (1|Clover), data=d6) # Correcting for iSize as fixed effect
  ycorr <- d6$growth_per_day - model.matrix( ~ InitialSize, data=d6) %*% fixef(fit) 
  d6$gpd_dryweight_Simplecor_asfixedeffect <- ycorr #this is the new corrected dry weight
  
}

traits_cor=cbind(d6$InitialSize,d6$growth_per_day,d6$gpd_dryweight_cor,d6$gpd_dryweight_Simplecor_asfixedeffect)
nrow(traits_cor) #replicates kept
colnames(traits_cor)=c("iSize","GPD_noCor","GPD_ResCor","GPD_FixedCor")


corrgram(traits_cor, order=F, panel=panel.pts, pch='.',
         upper.panel=panel.cor,cex=2)

plot(traits_cor[,1],traits_cor[,2],xlab="iSize", ylab="GPD", main="GPD no correction") #iSize and no corrected GPD
plot(traits_cor[,1],traits_cor[,3],xlab="iSize", ylab="GPD", main="GPD residual correction") #iSize and residual corrected GPD
plot(traits_cor[,1],traits_cor[,4], xlab="iSize", ylab="GPD", main="GPD fixed correction") #iSize and fixed corrected GPD

