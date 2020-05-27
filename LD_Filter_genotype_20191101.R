###############################################################################################
#           
#                 Script to replace a large number of SNPs with one consensus SNP pr. gene
#
###############################################################################################



library(stringr)
library(parallel)

# Read in data
{
setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GP_based_on_RNASNPs/20191101_Genotypefiltering/LD_filtering_based_on_genes/")
  
genefile=read.table("CHRPOS_unfiltered.genes.txt",sep="\t")
head(genefile)
dim(genefile) 
colnames(genefile)[1]="CHROM"
colnames(genefile)[4]="POS"

SNPswithgenes=subset(genefile,genefile[,19]==1)

genotype=read.table("MAFfilteredgenotypes23052019.csv",sep=",",header = T)
dim(genotype)

#Make a column in genefile with identifier for merging with genotypes
chridentifier=str_split_fixed(genefile$CHROM, "r", 2)[,2]
genefile$chridentifier =chridentifier
genefile$identifier=paste0(genefile$chridentifier,"-",genefile$POS)
genefile_small=genefile[,c(18,21)]
colnames(genefile_small)=c("gene","identifier")
genotype$identifier=paste0(genotype$CHROM,"-",genotype$POS)

All=merge(genotype,genefile_small,by.x = "identifier")
dim(All)

#Some SNPs are associated with more than one gene. 
length(unique(All$identifier)) #1,309,405 fits perfectly with the number of SNPs
nrow(All)/length(unique(All$identifier)) #In average every SNP has 1.6 genes associated with them

#Remove the outlier individuals!!

#outlier=which(colnames(All)=="Aalon_05")
#All=All[,-outlier]
#outlier=which(colnames(All)=="Banna_03")
#All=All[,-outlier]
#outlier=which(colnames(All)=="Ilona_07")
#All=All[,-outlier]
#outlier=which(colnames(All)=="Rling_02")
#All=All[,-outlier]
#outlier=which(colnames(All)=="Aoost_10")
#All=All[,-outlier]
#outlier=which(colnames(All)=="Volin_02")
#All=All[,-outlier]
}


# Prepare a gene column with chromosome followed by "-" and then the jg number identifying the gene
{
geneSHORT0=str_split_fixed(All$gene, ".jg", 2)[,2]
geneSHORT1=str_split_fixed(geneSHORT0, ".t", 2)[,1]
All$geneSHORT=geneSHORT1
All$geneSHORT=paste0(geneSHORT1,"-",All$CHROM)

All$geneSHORT[which(All$geneSHORT=="-0")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-1")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-2")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-3")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-4")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-5")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-6")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-7")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-8")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-9")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-10")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-11")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-12")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-13")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-14")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-15")]="No gene"
All$geneSHORT[which(All$geneSHORT=="-16")]="No gene"

length(unique(All$geneSHORT)) #40446 -1 = 40445 unique genes
}



#Not used 
#This part focus on keeping one SNP for each gene. 
# Actually one consensus SNP pr gene, lowering the number of SNPs to the number of genes plus the ones on chr0
{

#Calculates consensus
# The trick is now to choose one consensus SNP for all 78,906 different genes
# Make function that loops through each gene, takes the SNPs associated with it and calculate the consensus SNP 0,1,2 for each individual

#Consensus_for_a_gene=function(gene){
  #rows=which(All$geneSHORT==gene)
  #SNPs=All$identifier[rows]
  #NumberOfSNPs=length(SNPs)
  #Chr=All$CHROM[rows]
  #Pos=All$POS[rows]
  
  
  
  #for (i in seq(1:(NumberOfSNPs-1))){
    
    #if (NumberOfSNPs==1){
    #}
    
    #if (NumberOfSNPs>1){
      #if ((str_sub(SNPs[i],1,2)!=str_sub(SNPs[i+1],1,2))){
       # print("Problem. Gene not identified correct") 
      #}
      #}
  #}

  #CurrentSNPs=All[rows,4:129]
  #for (individual in seq(1:ncol(CurrentSNPs))){
    #zeros=length(which(CurrentSNPs[,individual]==0)) 
    #ones=length(which(CurrentSNPs[,individual]==1))
    #twos=length(which(CurrentSNPs[,individual]==2))
  
    #if(zeros+ones+twos!=NumberOfSNPs){
     # print("Problem: Genotypes found that were not 0,1 or 2")
    #}
    
    
    #if (which.max(c(zeros,ones,twos))==1){
      #most zero
     # print("most zero")
    #  Consensus=0
    #}
    
    #if (which.max(c(zeros,ones,twos))==2){
     # #most ones
      #print("most ones")
      #Consensus=1
  #  }
    
   # if (which.max(c(zeros,ones,twos))==3){
      #most ones
      #print("most twos")
      #Consensus=2
  #  }
  
   # if (c(zeros,ones,twos)[1]==c(zeros,ones,twos)[2]){
    #  #Equal numbers of zeros and ones
    #  print("equal numbers of ones and zeros")
    #  Consensus=0
  #  }
   # if (c(zeros,ones,twos)[1]==c(zeros,ones,twos)[3]){
      #Equal numbers of zeros and twos
    #  print("equal numbers of ones and twos")
     # Consensus=0
  #  }
   # if (c(zeros,ones,twos)[2]==c(zeros,ones,twos)[3]){
      #Equal numbers of ones and twos
    #  print("equal numbers of zeros and twos")
      #Consensus=2
    #}
    
    #Replace with consensus
    #CurrentSNPs[,individual]=Consensus
    
#  }
 # matrix=cbind(Chr[1],min(Pos),CurrentSNPs[1,],gene,NumberOfSNPs)
  #return(matrix)
#}

# Then apply the function to all gene categories
#{
#genelist=as.vector(unique(All$geneSHORT)) 
#results=mclapply(genelist[2:length(genelist)],Consensus_for_a_gene)
#results=Consensus_for_a_gene(All$geneSHORT[162609]) #testing, works

#}

# Make a genotype file with new consensus SNPs
#{
#x=rep(NA,length(results)*ncol(results[[1]]))
#genomatrix <- matrix(x,nrow=length(results),ncol=ncol(results[[1]]))
#colnames(genomatrix)=colnames(results[[1]])

#for (i in seq(1:length(results))){
 # genomatrix[i,]=unlist(results[[i]], use.names=FALSE)
#}

#head(genomatrix)

#}
  
  
#}


  
  
# Write a genotype file with genes for each SNP
{
Final=All[,c(2:135,137)]
#For some reasons some rows are duplicated. Remove those
deduped.data <- unique( Final[ , 1:ncol(Final) ] )
ordered=deduped.data[order(deduped.data$CHROM),]
write.table(ordered,"GenotypewithGenes_OutliersNOTremoved_20191108.csv",sep=",",col.names=TRUE,row.names=F,quote=F)


}

