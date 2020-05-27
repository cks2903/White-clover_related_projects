#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:59:46 2019

@author: CathrineKiel
"""

import pandas as pd
import numpy as np
from sys import argv
#  Marnis function
def calcGRM(SNP):
    """Code written by Marni"""
    N, M = SNP.shape #N is number of individuals, M is number of markers
    NORM = (SNP-np.mean(SNP, 0))/(np.std(SNP, 0)+0.000000000000001)
    return np.dot(NORM, NORM.T)/M


file=pd.read_csv(argv[1],header=0)
print("file read succesfully")

#delete the four accessions that probably belongs to another cultivar when looking at GRM
file = file.drop("Aalon_05", axis=1)
file = file.drop("Banna_03", axis=1)
file = file.drop("Ilona_07", axis=1)
file = file.drop("Rling_02", axis=1)
file = file.drop("Aoost_10", axis=1)
file = file.drop("Volin_02", axis=1)

row,col=file.shape

geno=np.array(file)
geno[0:1,]


# remove first two columns
geno1=np.copy(geno[:,2:col-2])
geno1[0:1,]
# Transpose
geno2=geno1.T
geno2[0:1,]

# replace all NA with mean genotype of that SNP
pandas_intermediate=pd.DataFrame(geno2)
col_mean=pandas_intermediate.mean(axis=0)
col_mean=np.array(col_mean)
inds=np.where(pd.isna(geno2))
geno2[inds] = np.take(col_mean, inds[1]) 
geno2=geno2.astype(float)
geno2[0:1,]

# Applying function, calculating GRM
GRM=calcGRM(geno2)
GRM[0:1,]


#Add header
names=file.columns[2:col-2]
names1=list(names)
GRM_df=pd.DataFrame(GRM,columns=names1)

GRM_df.to_csv(argv[2], index = None, header=True)
