##############################################################################
#             Checking heritability of growth periods in Clover. 
##############################################################################

# Load in data
{
  setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GP_based_on_RNASNPs/BacktoBasicYield_092019_GP/")
  library(lme4)
  library(BGLR)
  library("parallel")
  d=read.table("greenhouse_data_normalized_24092019.csv",header=T,sep=";")
  dim(d)
  
}

# A bit of filtering 
{
  before=dim(d)[1]
  d=d[-which(d$Rhizobium=="SM73"),]
  d=d[-which(d$Rhizobium=="NO"),]
  d=d[-which(d$n.stolons<4),]
  
  d$roundRep <- paste(d$Round, d$Replicate, sep='_')
  d=d[-which(d$roundRep=="1_2"),]
  
  after=dim(d)[1]
  before-after
  print(paste("removed",before-after,"observations,",after, "remains"),sep="")
}


# A function to return a matrix with all individuals and the growth pr. days during a certain period
growthprdaysinintervals<-function(start,end,data){
  matrix = matrix(,nrow = length(unique(data$Barcode)), ncol = 5)
  colnames(matrix)=c("Barcode","Start(days)","End(days)","number_of_days","GPD_in_interval")
  plants=unique(data$Barcode)
  for (plant in plants){
    name=paste(plant,"d",sep="")
    d2=subset(data,Barcode==plant)
    
    period_length=end-start+1

    days_matrix=rep(NA,100*2)
    days_matrix <- matrix(days_matrix,nrow=100,ncol=2)
    colnames(days_matrix)=c("day","index")
    for (day in seq(start,end,1)){
      dayidx=which(d2$NormTime==day+3)
      if (length(dayidx)!=0){
        days_matrix[day,]=c(day,dayidx)
      }
    }
    days_matrix=na.omit(days_matrix)
    
    NAS=period_length-nrow(days_matrix)    
    if (NAS<period_length*0.25){ #if more than 25% is missing in given growth interval, do not calculate GPD_interval
      start_idx=days_matrix[1,2]
      end_idx=days_matrix[nrow(days_matrix),2]
      reg_coef=lm(d2$Size[start_idx:end_idx]~(seq(start_idx:end_idx)))$coefficients[2]
      GPD_interval=reg_coef
      context=cbind(plant,start,end,period_length,GPD_interval)
      matrix[which(plants==plant),]=context
    }
    }
  return(matrix)
}


# Remove plants that do not grow from day 10 to 20
# These are considered dead or contaminated if they start growing after day 20
first20days=growthprdaysinintervals(10,20,d)
nrow(first20days)
hist(first20days[,5],breaks=200) 
mean(first20days[,5])
Dead=(which(first20days[,5]<100)) #These plants die
DeadBarcodes=first20days[Dead,1]
length(DeadBarcodes) #15 plants

# Take a look at the Barcodes being removed
x1821=subset(d,Barcode==1821) 
plot(x1821$NormTime,x1821$Size)

x1870=subset(d,Barcode==1870) 
plot(x1870$NormTime,x1870$Size)

x2277=subset(d,Barcode==2277) 
plot(x2277$NormTime,x2277$Size)

x3446=subset(d,Barcode==3446) 
plot(x3446$NormTime,x3446$Size)

x3502=subset(d,Barcode==3502) 
plot(x3502$NormTime,x3502$Size)

x3504=subset(d,Barcode==3504) 
plot(x3504$NormTime,x3504$Size)

x3522=subset(d,Barcode==3522) 
plot(x3522$NormTime,x3522$Size)

x3663=subset(d,Barcode==3663) 
plot(x3663$NormTime,x3663$Size)

x3869=subset(d,Barcode==3869) 
plot(x3869$NormTime,x3869$Size)

x3880=subset(d,Barcode==3880) 
plot(x3880$NormTime,x3880$Size)

x3979=subset(d,Barcode==3979) 
plot(x3979$NormTime,x3979$Size)

x4397=subset(d,Barcode==4397) 
plot(x4397$NormTime,x4397$Size)

x4462=subset(d,Barcode==4462) 
plot(x4462$NormTime,x4462$Size)

x4704=subset(d,Barcode==4704) 
plot(x4704$NormTime,x4704$Size)

x4902=subset(d,Barcode==4902) 
plot(x4902$NormTime,x4902$Size)


remove=which(d$Barcode %in% DeadBarcodes)
d1=d[-remove,]



### Remove plants than suddently drops in growth
plants=unique(d1$Barcode)
plants_that_behaves_weird=c()

for (i in 1:length(plants)){
  plant=plants[i]
  d2=subset(d1,Barcode==plant)
  
  for (j in 1:(nrow(d2)-1)){
    difference=d2$Size[j+1]-d2$Size[j]
    
    
    if (difference<(mean(d2$Size)*-0.9)){
      #print(plant)
      plants_that_behaves_weird[i]=plant
      break
    }
  }
  
}

plants_that_behaves_weird=na.omit(plants_that_behaves_weird)
length(plants_that_behaves_weird) #8 plants 
plants_that_behaves_weird=as.vector(plants_that_behaves_weird)
plants_that_behaves_weird

# Take a look at the Barcodes being removed
x3647=subset(d,Barcode==3647) 
plot(x3647$NormTime,x3647$Size)

x3650=subset(d,Barcode==3650) 
plot(x3650$NormTime,x3650$Size)

x4009=subset(d,Barcode==4009) 
plot(x4009$NormTime,x4009$Size)

x4209=subset(d,Barcode==4209) 
plot(x4209$NormTime,x4209$Size)

x5035=subset(d,Barcode==5035) 
plot(x5035$NormTime,x5035$Size)

x5132=subset(d,Barcode==5132) 
plot(x5132$NormTime,x5132$Size)

x5135=subset(d,Barcode==5135) 
plot(x5135$NormTime,x5135$Size)

x5168=subset(d,Barcode==5168) 
plot(x5168$NormTime,x5168$Size)


for (weirdos in plants_that_behaves_weird){
  indexes=which(d1$Barcode==weirdos)
  d1=d1[-indexes,]
}

length(unique(d1$Barcode)) #Still 2149 barcodes in dataset



# A function to calculate heritability for a certain growth period
HerCalc<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  
  # Now connect Barcodes with Clover, Rhizobium, EW, ES, RoundRep and calculate heritability for 1-5 days
  # Initial Size is an average of the size of the 5 first days
  # Correct for initial size
  
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=as.character(d1$Clover[which(d1$Barcode==dataframe$Barcode[i])[1]])
    dataframe$Rhizobium[i]=as.character(d1$Rhizobium[which(d1$Barcode==dataframe$Barcode[i])[1]])
    dataframe$EW[i]=d1$EW[which(d1$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d1$NS[which(d1$Barcode==dataframe$Barcode[i])[1]]
    dataframe$RoundRep[i]=d1$roundRep[which(d1$Barcode==dataframe$Barcode[i])[1]]
    #dataframe$InitialSize[i]=d1$InitialSize_avg5firstdays[which(d1$Barcode==dataframe$Barcode[i])[1]]
  }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + Rhizobium + factor(RoundRep) + (1|Clover), data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarClover=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarClover/(VarClover+VarResidual) 
  return(H2)
  
}










###############################################################################################
#
#           GPD in 5 days overlapping intervals and clover heritability calculation  
#
###############################################################################################

# Calculate GPD in 5 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix1=rep(NA,62*3)
matrix1 <- matrix(matrix1,nrow=62,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  Her=HerCalc(get(name))
  matrix1[i,]=cbind(c(start,end,Her))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


#write.table(matrix1,"H2_5days_intervals.csv",quote=F,row.names = F)


###############################################################################################
#
#           GPD in 10 days overlapping intervals and clover heritability calculation
#
###############################################################################################

# Calculate GPD in 10 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix2=rep(NA,62*3)
matrix2 <- matrix(matrix2,nrow=62,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=HerCalc(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix2,"H2_10days_intervals.csv",row.names=F,quote=F)

##############################################################################################
#
#           GPD in 15 days overlapping intervals and clover heritability calculation
#
###############################################################################################

# Calculate GPD in 15 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix3=rep(NA,62*3)
matrix3 <- matrix(matrix3,nrow=62,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=HerCalc(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])

plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


#write.table(matrix3,"H2_15days_intervals.csv",quote=F,row.names = F)

##############################################################################################
#
#           GPD in 20 days overlapping intervals and clover heritability calculation
#
###############################################################################################


# Calculate GPD in 20 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalc(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Clover heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


#write.table(matrix4,"H2_20days_intervals.csv",quote=F,row.names = F)

################################
# Get highly her phenotype     #
#day 7 to 26, median around day 16 or 17
head(day7to26)
phenotype=merge(d1,day7to26,by="Barcode")
length(unique(phenotype$Barcode))

phenotype1=phenotype[,c(1,5,7:12,15:17,22,26)]
nrow(phenotype1)

phenotype2=phenotype1[!duplicated(phenotype1$Barcode), ]       
nrow(phenotype2)

#write.table(phenotype2,"GPD_day7to26.csv",quote=F,row.names = F)


# Plot with all intervals
library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Clover heritability") +
  geom_hline(yintercept=0.45, linetype="dashed", color = "black", size=0.3) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))


# This plot is fine but I want to plot the median of an interval instead
# this will make it easier to compare the maximum of each curve, because they actually have maximums around the same time




# Calculate median for each period
matrix1_=data.frame(matrix1)
matrix1_=matrix1[1:30,]
matrix2_=data.frame(matrix2)
matrix2_=matrix2[1:30,]
matrix3_=data.frame(matrix3)
matrix3_=matrix3[1:30,]
matrix4_=data.frame(matrix4)
matrix4_=matrix4[1:30,]


matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
colnames(matrix_all)[2]="end.5"
colnames(matrix_all)[4]="end.10"
colnames(matrix_all)[6]="end.15"
colnames(matrix_all)[8]="end.20"
median.5=apply(matrix_all[,c("start","end.5")],1,median)
median.10=apply(matrix_all[,c("start","end.10")],1,median)
median.15=apply(matrix_all[,c("start","end.15")],1,median)
median.20=apply(matrix_all[,c("start","end.20")],1,median)


newmatrix5days=cbind(median.5,matrix_all$H2_5days)
colnames(newmatrix5days)=c("median","H2_5days")
newmatrix10days=cbind(median.10,matrix_all$H2_10days)
colnames(newmatrix10days)=c("median","H2_10days")
newmatrix15days=cbind(median.15,matrix_all$H2_15days)
colnames(newmatrix15days)=c("median","H2_15days")
newmatrix20days=cbind(median.20,matrix_all$H2_20days)
colnames(newmatrix20days)=c("median","H2_20days")


matrix_all_=rbind(newmatrix5days,newmatrix10days,newmatrix15days,newmatrix20days)
matrix_all_=as.data.frame(matrix_all_)
matrix_all_$group[1:30]="5 days intervals"
matrix_all_$group[31:60]="10 days intervals"
matrix_all_$group[61:90]="15 days intervals"
matrix_all_$group[91:nrow(matrix_all_)]="20 days intervals"

matrix_all_=matrix_all_[order(matrix_all_[,1]),] 

#you can remove all median values above 30
matrix_all_=matrix_all_[-which(matrix_all_$median>30),]


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Clover heritability") +
  geom_hline(yintercept=0.45, linetype="dashed", color = "black", size=0.3) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))


# Now try to add error bars to the point
matrix_all_$xmin=NA
matrix_all_$xmax=NA



growthinterval1=which(matrix_all_$group=="5 days intervals")
matrix_all_$xmin[growthinterval1]=matrix_all_$median[growthinterval1]-2
matrix_all_$xmax[growthinterval1]=matrix_all_$median[growthinterval1]+2

growthinterval2=which(matrix_all_$group=="10 days intervals")
matrix_all_$xmin[growthinterval2]=matrix_all_$median[growthinterval2]-4.5
matrix_all_$xmax[growthinterval2]=matrix_all_$median[growthinterval2]+4.5

growthinterval3=which(matrix_all_$group=="15 days intervals")
matrix_all_$xmin[growthinterval3]=matrix_all_$median[growthinterval3]-7
matrix_all_$xmax[growthinterval3]=matrix_all_$median[growthinterval3]+7

growthinterval4=which(matrix_all_$group=="20 days intervals")
matrix_all_$xmin[growthinterval4]=matrix_all_$median[growthinterval4]-9.5
matrix_all_$xmax[growthinterval4]=matrix_all_$median[growthinterval4]+9.5


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  geom_errorbarh(data=matrix_all_,aes(xmin = xmin,xmax = xmax),color="gray",alpha=0.4,height=0, size=0.5) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Clover heritability") +
  geom_hline(yintercept=0.45, linetype="dashed", color = "black", size=0.3) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold")) 


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  geom_errorbarh(data=matrix_all_,aes(color=group,alpha=0.01,xmin = xmin,xmax = xmax),height=0, size=0.3) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Clover heritability") +
  geom_hline(yintercept=0.45, linetype="dashed", color = "black", size=0.3) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Clover heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold")) 













########################################################################################################################################
#                                                                                                                                      #
#                   Now calculate rhizobium heritability instead of clover                                                             #
#                                                                                                                                      #
########################################################################################################################################
{
#Filtering as in Clover heritability calculations


HerCalcRhiz<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  
  # Now connect Barcodes with Clover, Rhizobium, EW, ES, RoundRep and calculate heritability for 1-5 days
  # Initial Size is an average of the size of the 5 first days
  # Correct for initial size
  
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=as.character(d1$Clover[which(d1$Barcode==dataframe$Barcode[i])[1]])
    dataframe$Rhizobium[i]=as.character(d1$Rhizobium[which(d1$Barcode==dataframe$Barcode[i])[1]])
    dataframe$EW[i]=d1$EW[which(d1$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d1$NS[which(d1$Barcode==dataframe$Barcode[i])[1]]
    dataframe$RoundRep[i]=d1$roundRep[which(d1$Barcode==dataframe$Barcode[i])[1]]
    #dataframe$InitialSize[i]=d1$InitialSize_avg5firstdays[which(d1$Barcode==dataframe$Barcode[i])[1]]
  }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + (1|Rhizobium) + factor(RoundRep) + Clover, data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarClover=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarClover/(VarClover+VarResidual) 
  return(H2)
  
}





# Calculate GPD in 5 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix1=rep(NA,62*3)
matrix1 <- matrix(matrix1,nrow=62,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  Her=HerCalcRhiz(get(name))
  matrix1[i,]=cbind(c(start,end,Her))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix1,"Rhiz_H2_5days_intervals.csv",quote=F,row.names = F)



# Calculate GPD in 10 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix2=rep(NA,62*3)
matrix2 <- matrix(matrix2,nrow=62,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=HerCalcRhiz(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix2,"Rhiz_H2_10days_intervals.csv",row.names=F,quote=F)



# Calculate GPD in 15 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix3=rep(NA,62*3)
matrix3 <- matrix(matrix3,nrow=62,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=HerCalcRhiz(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])

plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


#write.table(matrix3,"Rhiz_H2_15days_intervals.csv",quote=F,row.names = F)





# Calculate GPD in 20 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1))
}


# Heritability calculation; applying the function to different growth intervals

matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcRhiz(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"Rhiz_H2_20days_intervals.csv",quote=F,row.names = F)

# Plot with all intervals
library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))

# Calculate median for each period
matrix1_=data.frame(matrix1)
matrix1_=matrix1[1:30,]
matrix2_=data.frame(matrix2)
matrix2_=matrix2[1:30,]
matrix3_=data.frame(matrix3)
matrix3_=matrix3[1:30,]
matrix4_=data.frame(matrix4)
matrix4_=matrix4[1:30,]


matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
colnames(matrix_all)[2]="end.5"
colnames(matrix_all)[4]="end.10"
colnames(matrix_all)[6]="end.15"
colnames(matrix_all)[8]="end.20"
median.5=apply(matrix_all[,c("start","end.5")],1,median)
median.10=apply(matrix_all[,c("start","end.10")],1,median)
median.15=apply(matrix_all[,c("start","end.15")],1,median)
median.20=apply(matrix_all[,c("start","end.20")],1,median)


newmatrix5days=cbind(median.5,matrix_all$H2_5days)
colnames(newmatrix5days)=c("median","H2_5days")
newmatrix10days=cbind(median.10,matrix_all$H2_10days)
colnames(newmatrix10days)=c("median","H2_10days")
newmatrix15days=cbind(median.15,matrix_all$H2_15days)
colnames(newmatrix15days)=c("median","H2_15days")
newmatrix20days=cbind(median.20,matrix_all$H2_20days)
colnames(newmatrix20days)=c("median","H2_20days")


matrix_all_=rbind(newmatrix5days,newmatrix10days,newmatrix15days,newmatrix20days)
matrix_all_=as.data.frame(matrix_all_)
matrix_all_$group[1:30]="5 days intervals"
matrix_all_$group[31:60]="10 days intervals"
matrix_all_$group[61:90]="15 days intervals"
matrix_all_$group[91:nrow(matrix_all_)]="20 days intervals"

matrix_all_=matrix_all_[order(matrix_all_[,1]),] 

#you can remove all median values above 30
matrix_all_=matrix_all_[-which(matrix_all_$median>30),]


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))


# Now try to add error bars to the point
matrix_all_$xmin=NA
matrix_all_$xmax=NA



growthinterval1=which(matrix_all_$group=="5 days intervals")
matrix_all_$xmin[growthinterval1]=matrix_all_$median[growthinterval1]-2
matrix_all_$xmax[growthinterval1]=matrix_all_$median[growthinterval1]+2

growthinterval2=which(matrix_all_$group=="10 days intervals")
matrix_all_$xmin[growthinterval2]=matrix_all_$median[growthinterval2]-4.5
matrix_all_$xmax[growthinterval2]=matrix_all_$median[growthinterval2]+4.5

growthinterval3=which(matrix_all_$group=="15 days intervals")
matrix_all_$xmin[growthinterval3]=matrix_all_$median[growthinterval3]-7
matrix_all_$xmax[growthinterval3]=matrix_all_$median[growthinterval3]+7

growthinterval4=which(matrix_all_$group=="20 days intervals")
matrix_all_$xmin[growthinterval4]=matrix_all_$median[growthinterval4]-9.5
matrix_all_$xmax[growthinterval4]=matrix_all_$median[growthinterval4]+9.5


ggplot(data=matrix_all_,aes(x=median, y=H2_5days,group=group,fill=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group)) +
  geom_errorbarh(data=matrix_all_,aes(xmin = xmin,xmax = xmax),color="gray",alpha=0.4,height=0, size=0.5) +
  scale_color_manual(values=c("#1B9E77","#386CB0","#FB8072","#984EA3")) +
  xlab("Middle of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold")) 


}

#############################################################################################
### Now calculate Rhizobium Heritability by grouping Rhizobium genotypes into genospecies####
#############################################################################################

#Merge with information about genospecies
{
genospe=read.table("genospecies.txt",sep="\t",head=T)
colnames(genospe)[1]="Rhizobium"

d1_genospe=merge(d1,genospe,by="Rhizobium") #Full table with genospecies information
}

{
HerCalcGenospecies<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=as.character(d1_genospe$Clover[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$Genospecies[i]=as.character(d1_genospe$Genospecies[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$EW[i]=d1_genospe$EW[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d1_genospe$NS[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$RoundRep[i]=d1_genospe$roundRep[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
  }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + (1|Genospecies) + factor(RoundRep) + Clover, data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarClover=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarClover/(VarClover+VarResidual) 
  return(H2)
  
}


# Heritability calculation; applying the function to different growth intervals

matrix1=rep(NA,62*3)
matrix1 <- matrix(matrix1,nrow=62,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  Her=HerCalcGenospecies(get(name))
  matrix1[i,]=cbind(c(start,end,Her))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix1,"Genospecies_H2_5days_intervals.csv",quote=F,row.names = F)


# Heritability calculation; applying the function to different growth intervals

matrix2=rep(NA,62*3)
matrix2 <- matrix(matrix2,nrow=62,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=HerCalcGenospecies(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix2,"Genospecies_H2_10days_intervals.csv",row.names=F,quote=F)



matrix3=rep(NA,62*3)
matrix3 <- matrix(matrix3,nrow=62,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=HerCalcGenospecies(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])

plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


#write.table(matrix3,"Genospecies_H2_15days_intervals.csv",quote=F,row.names = F)


matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcGenospecies(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"Genospecies_H2_20days_intervals.csv",quote=F,row.names = F)

# Plot with all intervals
library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))


}

##################################################################################################
### Calculate Rhizobium Heritability by grouping Rhizobium genotypes into Nodphyl genospecies ####
##################################################################################################

{

HerCalcNodphyl<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=as.character(d1_genospe$Clover[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$Nodphyl_fixT[i]=as.character(d1_genospe$Nodphyl_fixT[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$EW[i]=d1_genospe$EW[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d1_genospe$NS[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$RoundRep[i]=d1_genospe$roundRep[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
  }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + (1|Nodphyl_fixT) + factor(RoundRep) + Clover, data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarClover=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarClover/(VarClover+VarResidual) 
  return(H2)
  
}


# Heritability calculation; applying the function to different growth intervals

matrix1=rep(NA,62*3)
matrix1 <- matrix(matrix1,nrow=62,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  Her=HerCalcNodphyl(get(name))
  matrix1[i,]=cbind(c(start,end,Her))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix1,"Nodphyl_Rhiz_H2_5days_intervals.csv",quote=F,row.names = F)


# Heritability calculation; applying the function to different growth intervals

matrix2=rep(NA,62*3)
matrix2 <- matrix(matrix2,nrow=62,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=HerCalcNodphyl(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix2,"Nodphyl_Rhiz_H2_10days_intervals.csv",row.names=F,quote=F)



matrix3=rep(NA,62*3)
matrix3 <- matrix(matrix3,nrow=62,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=HerCalcNodphyl(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])

plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))


#write.table(matrix3,"Nodphyl_Rhiz_H2_15days_intervals.csv",quote=F,row.names = F)


matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcNodphyl(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"Nodphyl_Rhiz_H2_20days_intervals.csv",quote=F,row.names = F)

# Plot with all intervals
library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability (Nodphyl) for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))

}

##################################################################################################
### Calculate Rhizobium Heritability using only Aearl_0749                                    ####
##################################################################################################

#Aearl_0749 has been inoculated with all rhizobium strains meaning that heritability of rhizobium
# can be calculated here without adding clover as a random effect

{
Aearl_0749=d1_genospe[which(d1_genospe$Clover=="Aearl_0749"),]
head(Aearl_0749)
nrow(Aearl_0749) #21009



# Calculate GPD in 5 days overlapping intervals 
for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,Aearl_0749))
}

length(unique(day1to5[,1])) #385, 163 unique rhizobium strains


HerCalcRhizAearl0749<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  
  # Now connect Barcodes with Clover, Rhizobium, EW, ES, RoundRep and calculate heritability for 1-5 days
  # Initial Size is an average of the size of the 5 first days
  # Correct for initial size
  
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=as.character(d1_genospe$Clover[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$Rhizobium[i]=as.character(d1_genospe$Rhizobium[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$EW[i]=d1_genospe$EW[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d1_genospe$NS[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$RoundRep[i]=d1_genospe$roundRep[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    #dataframe$InitialSize[i]=d1$InitialSize_avg5firstdays[which(d1$Barcode==dataframe$Barcode[i])[1]]
  }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + (1|Rhizobium) + factor(RoundRep), data=dataframe)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarClover=re_dat[1,'vcov']
  VarResidual=re_dat[2,'vcov']
  H2=VarClover/(VarClover+VarResidual) 
  return(H2)
  
}

matrix1=rep(NA,62*3)
matrix1 <- matrix(matrix1,nrow=62,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  Her=HerCalcRhizAearl0749(get(name))
  matrix1[i,]=cbind(c(start,end,Her))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix1,"Aearl_0749_rhiz_H2_5days_intervals.csv",quote=F,row.names = F)



matrix2=rep(NA,62*3)
matrix2 <- matrix(matrix2,nrow=62,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=HerCalcRhizAearl0749(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix2,"Aearl_0749_rhiz_H2_10days_intervals.csv",quote=F,row.names = F)


matrix3=rep(NA,62*3)
matrix3 <- matrix(matrix3,nrow=62,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=HerCalcRhizAearl0749(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])

plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix3,"Aearl_0749_rhiz_H2_15days_intervals.csv",quote=F,row.names = F)



matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcRhizAearl0749(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"Aearl_0749_rhiz_H2_20days_intervals.csv",quote=F,row.names = F)


library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))

}

##################################################################################################
### Calculate Rhizobium  (genospecies) Heritability using only Aearl_0749                     ####
##################################################################################################

{
#Aearl_0749 has been inoculated with all rhizobium strains meaning that heritability of rhizobium
# can be calculated here without adding clover as a random effect

  HerCalcGenospecies_Aearl0749<-function(dataframe){
    dataframe=as.data.frame(dataframe)
    for (i in 1:nrow(dataframe)){
      dataframe$Clover[i]=as.character(Aearl_0749$Clover[which(Aearl_0749$Barcode==dataframe$Barcode[i])[1]])
      dataframe$Genospecies[i]=as.character(Aearl_0749$Genospecies[which(Aearl_0749$Barcode==dataframe$Barcode[i])[1]])
      dataframe$EW[i]=Aearl_0749$EW[which(Aearl_0749$Barcode==dataframe$Barcode[i])[1]]
      dataframe$NS[i]=Aearl_0749$NS[which(Aearl_0749$Barcode==dataframe$Barcode[i])[1]]
      dataframe$RoundRep[i]=Aearl_0749$roundRep[which(Aearl_0749$Barcode==dataframe$Barcode[i])[1]]
    }
    dataframe=na.omit(dataframe)
    # Check heritabilities for all growth periods
    Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + (1|Genospecies) + factor(RoundRep), data=dataframe)
    summary(Correctedforallfixed)
    re_dat = as.data.frame(VarCorr(Correctedforallfixed))
    VarClover=re_dat[1,'vcov']
    VarResidual=re_dat[2,'vcov']
    H2=VarClover/(VarClover+VarResidual) 
    return(H2)
    
  }


matrix1=rep(NA,62*3)
matrix1 <- matrix(matrix1,nrow=62,ncol=3)
colnames(matrix1)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+4
  name=paste("day",i,"to",i+4,sep="")
  get(name)
  Her=HerCalcGenospecies_Aearl0749(get(name))
  matrix1[i,]=cbind(c(start,end,Her))
  
}

matrix1
matrix1=na.omit(matrix1)
plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1))
cor(matrix1[,1],matrix1[,3])

plot(matrix1[,1],matrix1[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 5 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix1,"Aearl_0749_rhizgenospecies_H2_5days_intervals.csv",quote=F,row.names = F)

matrix2=rep(NA,62*3)
matrix2 <- matrix(matrix2,nrow=62,ncol=3)
colnames(matrix2)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+9
  name=paste("day",i,"to",i+9,sep="")
  get(name)
  Her=HerCalcGenospecies_Aearl0749(get(name))
  matrix2[i,]=cbind(c(start,end,Her))
  
}

matrix2
matrix2=na.omit(matrix2)
plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1))
cor(matrix2[,1],matrix2[,3])

plot(matrix2[,1],matrix2[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 10 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix2,"Aearl_0749_rhizgenospecies_H2_10days_intervals.csv",quote=F,row.names = F)


matrix3=rep(NA,62*3)
matrix3 <- matrix(matrix3,nrow=62,ncol=3)
colnames(matrix3)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",i,"to",i+14,sep="")
  get(name)
  Her=HerCalcGenospecies_Aearl0749(get(name))
  matrix3[i,]=cbind(c(start,end,Her))
  
}

matrix3
matrix3=na.omit(matrix3)
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1))
cor(matrix3[,1],matrix3[,3])
plot(matrix3[,1],matrix3[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 15 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix3,"Aearl_0749_rhizgenospecies_H2_15days_intervals.csv",quote=F,row.names = F)



matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcGenospecies_Aearl0749(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium heritability (genospecies) of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"Aearl_0749_rhizgenospecies_H2_20days_intervals.csv",quote=F,row.names = F)


library(ggplot2)
library(reshape2)

matrix1_=data.frame(matrix1)
matrix2_=data.frame(matrix2)
matrix3_=data.frame(matrix3)
matrix4_=data.frame(matrix4)

matrix_all=merge(matrix1_,matrix2_,by= "start")
matrix_all=merge(matrix_all,matrix3_,by="start")
matrix_all=merge(matrix_all,matrix4_,by="start")
colnames(matrix_all)[3]="H2_5days"
colnames(matrix_all)[5]="H2_10days"
colnames(matrix_all)[7]="H2_15days"
colnames(matrix_all)[9]="H2_20days"
matrix_all[,2]=NULL
matrix_all[,3]=NULL
matrix_all[,4]=NULL
matrix_all[,5]=NULL
head(matrix_all)

premelted=matrix_all[1:30,]
melted=melt(premelted,id.var="start")


ggplot(data=melted,aes(x=start, y=value,group=variable,fill=variable)) +
  geom_line(aes(color=variable))+
  geom_point(aes(color=variable)) +
  scale_color_manual(values=c("#984EA3","#1B9E77","#386CB0","#FB8072")) +
  xlab("Start of interval") + ylab("Rhizobium heritability") +
  geom_hline(yintercept=0.20, linetype="dashed", color = "white", size=0.3,alpha=0) +
  scale_y_continuous(breaks = seq(from = 0, to = 1.00, by = 0.05)) +
  scale_x_continuous(breaks = seq(from = 0, to = 30, by = 2)) +
  ggtitle("Rhizobium heritability Aearl_0749 for different growth intervals") +
  theme( panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "grey50",size=0.05),
         plot.title=element_text(color="black", size=16, face="bold"))

}


########################################################################################################################################
#                                                                                                                                      #
#                   Clover -Rhizobium interaction term                                                                                 #
#                                                                                                                                      #
########################################################################################################################################

d1_genospe$Clo_Rhi_int = paste(d1_genospe$Clover, d1_genospe$Rhizobium, sep=':')


HerCalcInteraction<-function(dataframe){
  dataframe=as.data.frame(dataframe)
  
  # Now connect Barcodes with Clover, Rhizobium, EW, ES, RoundRep and calculate heritability for 1-5 days
  # Initial Size is an average of the size of the 5 first days
  # Correct for initial size
  
  for (i in 1:nrow(dataframe)){
    dataframe$Clover[i]=as.character(d1_genospe$Clover[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$Rhizobium[i]=as.character(d1_genospe$Rhizobium[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    dataframe$EW[i]=d1$EW[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$NS[i]=d1$NS[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$RoundRep[i]=d1_genospe$roundRep[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]]
    dataframe$Clo_Rhi_int[i]=as.character(d1_genospe$Clo_Rhi_int[which(d1_genospe$Barcode==dataframe$Barcode[i])[1]])
    }
  dataframe=na.omit(dataframe)
  # Check heritabilities for all growth periods
  Correctedforallfixed <- lmer(GPD_in_interval ~ factor(EW) + factor(NS) + (1|Clo_Rhi_int) + factor(RoundRep) + (1|Clover), data=dataframe)
  anova(Correctedforallfixed)
  summary(Correctedforallfixed)
  re_dat = as.data.frame(VarCorr(Correctedforallfixed))
  VarInt=re_dat[1,'vcov']
  VarClover=re_dat[2,'vcov']
  VarResidual=re_dat[3,'vcov']
  H2=VarInt/(VarClover+VarInt+VarResidual) 
  return(H2)
}

# Calculate GPD in 15 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+14
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1_genospe))
}

# Heritability calculation; applying the function to different growth intervals

matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcInteraction(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium:Clover interaction heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium:Clover heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"RhizCloverInt_H2_20days_intervals.csv",quote=F,row.names = F)




# Calculate GPD in 20 days overlapping intervals 

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",start,"to",end,sep="")
  assign(name,growthprdaysinintervals(start,end,d1_genospe))
}

# Heritability calculation; applying the function to different growth intervals

matrix4=rep(NA,62*3)
matrix4 <- matrix(matrix4,nrow=62,ncol=3)
colnames(matrix4)=c("start","end","H2")

for (i in 1:62){
  start=i
  end=i+19
  name=paste("day",i,"to",i+19,sep="")
  get(name)
  Her=HerCalcInteraction(get(name))
  matrix4[i,]=cbind(c(start,end,Her))
  
}

matrix4
matrix4=na.omit(matrix4)
plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium:Clover interaction heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1))
cor(matrix4[,1],matrix4[,3])

plot(matrix4[,1],matrix4[,3],xlab="Start day",ylab="Heritability",main="Rhizobium:Clover heritability of GPD for growth periods of 20 days",type="l",col="blue",ylim=c(0,1),xlim=c(0,30))

#write.table(matrix4,"RhizCloverInt_H2_20days_intervals.csv",quote=F,row.names = F)
