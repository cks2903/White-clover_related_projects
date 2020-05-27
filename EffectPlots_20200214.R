
# Get genotypes
geno=read.table("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GP_based_on_RNASNPs/LD_filtering_of_RNAseq/New_LD_filtering_within_genes_20191125/Geno_LDfilteredwithingenes_20191211.csv",sep=",",header=T)

setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GWAS results/GWAS_20200206/")

dim(geno)

pheno=read.table("Phenotypes_GPD.csv",header=T,sep=";")
dim(pheno)

pheno_int=read.table("GPDint_day7to26.csv",header=T,sep=",")
dim(pheno_int)

geno1=geno[,which(colnames(geno) %in% pheno$Accession==T)]
ncol(geno1)

geno2=cbind(geno[,1:2],geno1)
ncol(geno2)

# effect plots

# Top snp for day7to26 gpd
snp1=geno2[which(geno2$POS=="6854555"),]
snp1_=snp1[3:length(snp1)]
snp1_=t(snp1_)
snp1_=data.matrix(snp1_)
snp1_=as.data.frame(snp1_)
library(data.table)
setDT(snp1_, keep.rownames = TRUE)[]
colnames(snp1_)=c("Accession","genotype")


plot1=merge(snp1_,pheno,by="Accession")
nrow(plot1)
boxplot(plot1$GPD_Fixcor_iSize_Full~plot1$genotype,main="Chr1_Pos311818694",xlab="Genotype",ylab="GPD_Day7to26")
stripchart(GPD_Fixcor_iSize_Full ~ genotype, vertical = TRUE, data = plot1, 
                                                    method = "jitter", add = TRUE, pch = 20, col = 'orange') #save as 600 x 550 png

# Top snp for GPD_FixedCor
snp2=geno2[which(geno2$POS=="11818694"),]
snp2_=snp2[3:length(snp2)]
snp2_=t(snp2_)
snp2_=data.matrix(snp2_)
snp2_=as.data.frame(snp2_)
setDT(snp2_, keep.rownames = TRUE)[]
colnames(snp2_)=c("Accession","genotype")
plot2=cbind(snp2_,pheno_int)
nrow(plot2)


boxplot(plot2$pheno_int~plot2$genotype,main="Chr1_Pos11818694",xlab="Genotype",ylab="GPD_FixedCor")
stripchart(GPD_Fixcor_iSize_Full ~ genotype, vertical = TRUE, data = plot2, 
           method = "jitter", add = TRUE, pch = 20, col = 'orange') #save as 600 x 550 png

snp3=geno2[which(geno2$POS=="41740701"),]
snp3_=snp3[3:length(snp3)]
snp3_=t(snp3_)
snp3_=data.matrix(snp3_)
snp3_=as.data.frame(snp3_)
setDT(snp3_, keep.rownames = TRUE)[]
colnames(snp3_)=c("Accession","genotype")
plot3=cbind(snp3_,pheno)
nrow(plot3)


boxplot(plot3$GPD_Fixcor_iSize_Full~plot3$genotype,main="Chr16_Pos41740701",xlab="Genotype",ylab="GPD_FixedCor")
stripchart(GPD_Fixcor_iSize_Full ~ genotype, vertical = TRUE, data = plot3, 
           method = "jitter", add = TRUE, pch = 20, col = 'orange') #save as 600 x 550 png

snp4=geno2[which(geno2$POS=="75462622"),]
snp4_=snp4[3:length(snp4)]
snp4_=t(snp4_)
snp4_=data.matrix(snp4_)
snp4_=as.data.frame(snp4_)
setDT(snp4_, keep.rownames = TRUE)[]
colnames(snp4_)=c("Accession","genotype")
plot4=cbind(snp4_,pheno)
nrow(plot4)


boxplot(plot4$GPD_Fixcor_iSize_Full~plot4$genotype,main="Chr1_Pos75462622",xlab="Genotype",ylab="GPD_FixedCor")
stripchart(GPD_Fixcor_iSize_Full ~ genotype, vertical = TRUE, data = plot4, 
           method = "jitter", add = TRUE, pch = 20, col = 'orange') #save as 600 x 550 png



