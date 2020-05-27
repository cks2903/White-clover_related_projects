########################################################################################################
# MANHATTAN PLOTS GPD corrected data of 112 individuals (6 outliers removed) using LD-filtered SNPs    #
########################################################################################################

#source("https://bioconductor.org/biocLite.R")
#biocLite("IHW")
library("IHW")

setwd("/Volumes/ST_MBG-PMg/Cathrine/Nchain/Files_for_GWAS_and_genomic_prediction/GWAS results/GWAS_20200206/Results")
d=read.table("_pid4_GPD_ResiCor_iSize_Full_emmax_none_t75.pvals",head=T,sep="\t")

dim(d) #984,275 SNPs

library(qqman)

zero=which(d$chromosomes==0)
d1=d[-zero,]

MAFfilter=which(d1$macs<6)
d2=d1[-MAFfilter,]

Bonferoni=0.05/nrow(d2)
line=-log10(Bonferoni)
pvalues=as.vector(as.numeric(d$scores))
#BH threshold: Find the largest p-value smaller than the critical value
pvalues_ordered=sort(pvalues, decreasing = FALSE)
Criticalvalues=(1/length(pvalues_ordered)*0.05)

pvalues_ordered[1]<Criticalvalues #None of the SNPs are significant 

BenjaminiHochberg=get_bh_threshold(pvalues, 0.1, mtests = length(pvalues))


colnames(d2)=c("CHR","BP","P","MAFs","MACs","genotype_var_perc")
manhattan(d2,ylim=c(0,9))  #save as 2000 x 500 png

manhattan(subset(d2, CHR == 2), main = "Chr 2") # red line is Bonferoni threshold

manhattan(subset(d2, CHR == 7), main = "Chr 7") #Save as 1000x500 png

manhattan(subset(d2, CHR == 1), main = "Chr 1") #Save as 1000x500 png









