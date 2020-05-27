#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:59:46 2019

@author: CathrineKiel
"""

import os
import pandas as pd
import numpy as np
from sys import argv
#os.chdir("/Users/CathrineKiel/Desktop/GWAS_TRAINING_NEW_SNPS")

file=pd.read_csv(argv[1],header=0)
#"_pid1_Yellow.tip_emmax_none_t75_round2.pvals"
# first filter out MAF<0.05

row,col=file.shape

remove=[]
for i in range(0,row):
    
    if file.mafs[i]<0.05:
        remove.append(i)

len(remove)
MAFfiltered=file.drop(file.index[[remove]])

len(MAFfiltered)

# Print a list of the 200 SNPs to make GRM from
MAFfiltered.sort_values('scores',inplace=True)
top200=MAFfiltered.head(200)
headers=list(top200.columns)
headers=str(headers)

np.savetxt(argv[3],top200,fmt='%s',delimiter=',',header=headers)
#YellowTipRound2top200snps.txt"

# Now load in genotype file and remove everything that is not the top 200 SNPs
# We need to have the genotype file with all individuals, not just training
geno = pd.read_csv(argv[2],header=0)
#Cloverimputed2019final.csv

#If gpd remember to delete the individuals with no phenotype records and the 6 outlier individuals
dropcolumns=["Banna_05","Ccyma_05","Ccyma_07","Ccyma_08","Ctain_08","Kdike_05","Llanc_01","Llanc_02","Llanc_10","Mrida_03","Volin_06","Volin_08" , "Volin_10","Aalon_05","Banna_03","Ilona_07","Rling_02","Aoost_10","Volin_02","Aaran_02"]
geno=geno.drop(dropcolumns,axis=1)


looking_for=top200.ix[:,0:2]

geno_quick=np.array(geno) #converted genotype dataframe to numpy array for quicker looping.

indexes_needed_in_geno=[]

for i in range(0,200):
    current_suspect=looking_for.loc[looking_for.index[i],'positions']
    index=np.where(geno_quick == current_suspect)
    
    if len(index[0])==0:
        print(current_suspect,"not in genotype file")
    
    if len(index[0])==1:
        indexes_needed_in_geno.append(index[0][0])
        print(current_suspect,"found")
    
    if len(index[0])>1: #this is if the position ID is not unique
        chromosome=looking_for.loc[looking_for.index[i],'chromosomes']
        
        lst1=list(np.where(geno_quick == current_suspect))
        lst2=list(np.where(geno_quick == chromosome))
        
        c = [x for x in lst1[0] if x in lst2[0]]
        indexes_needed_in_geno.append(c[0])


indexes_needed_in_geno
len(indexes_needed_in_geno)
newgeno = geno_quick[indexes_needed_in_geno]
newgeno.shape
newgeno=np.array(newgeno)

# Then make GRM based on top 200, using Marnis function
def calcGRM(SNP):
    """Code written by Marni"""
    N, M = SNP.shape #N is number of individuals, M is number of markers
    NORM = (SNP-np.mean(SNP, 0))/(np.std(SNP, 0)+0.000000000000001)
    return np.dot(NORM, NORM.T)/M

# remove first two columns
geno1=np.copy(newgeno[:,2:(newgeno.shape[1]-1)])

# Transpose
geno2=geno1.T
geno2=geno2.astype(float)

# Applying function, calculating GRM
GRM=calcGRM(geno2)

#Add header
names=geno.columns[2:(newgeno.shape[1]-1)]
names1=list(names)
GRM_df=pd.DataFrame(GRM,columns=names1)

# Save file
# Write GRM based on top 200 significant GWAS SNPs in trn pop
GRM_df.to_csv(argv[4], index = None, header=True)
