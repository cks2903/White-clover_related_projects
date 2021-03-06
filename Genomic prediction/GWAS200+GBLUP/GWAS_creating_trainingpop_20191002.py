#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A script to take testing population and remove those from the genotype file
"""

#import os
#os.chdir("/Users/CathrineKiel/Desktop/GWAS_GP_20191002")


import sys
import pandas as pd

def read_as_list_of_lists(filename):
    _file = open(filename)
    data = []
    temp = []
    for line in _file:
        line = line.strip()
        if line=="": continue
        if line[:2]=="[[":
            if len(temp)!=0:
                data.append(temp)
                temp = []
            continue
        if line[0]=="[":
            for item in line.split(" ")[1:]:
                temp.append(item.replace('"', ''))
    data.append(temp)

    return data


if __name__=="__main__":
    testing_pop = read_as_list_of_lists(sys.argv[1])
    print(testing_pop)

 
group1tst=testing_pop[0]
group2tst=testing_pop[1]
group3tst=testing_pop[2]
group4tst=testing_pop[3]
group5tst=testing_pop[4]
group6tst=testing_pop[5]






file=pd.read_csv(sys.argv[2],header=0)
file.columns
drop=["identifier"]
file=file.drop(drop,axis=1)

#remove individuals which we don't have any phenotype data on and outliers
dropcolumns=["Aalon_05","Aoost_10","Banna_03","Banna_05","Ccyma_05","Ccyma_07","Ccyma_08","Ctain_08","Ilona_07","Kdike_05","Llanc_01","Llanc_02","Llanc_10","Mrida_03","Rling_02","Volin_02","Volin_06","Volin_08" , "Volin_10"]
file=file.drop(dropcolumns,axis=1)
file.shape

# Remove testing pop1
traininggeno1=file.drop(group1tst, axis=1)
traininggeno1.columns

# Remove testing pop2
traininggeno2=file.drop(group2tst, axis=1)
traininggeno2.columns

# Remove testing pop3
traininggeno3=file.drop(group3tst, axis=1)
traininggeno3.columns

# Remove testing pop4
traininggeno4=file.drop(group4tst, axis=1)
traininggeno4.columns

# Remove testing pop5
traininggeno5=file.drop(group5tst, axis=1)
traininggeno5.columns

# Remove testing pop6
traininggeno6=file.drop(group6tst, axis=1)
traininggeno6.columns

traininggeno1.to_csv("GWAStraining1_GPD.csv", index = None, header= 1)
traininggeno2.to_csv("GWAStraining2_GPD.csv", index = None, header= 1)
traininggeno3.to_csv("GWAStraining3_GPD.csv", index = None, header= 1)
traininggeno4.to_csv("GWAStraining4_GPD.csv", index = None, header= 1)
traininggeno5.to_csv("GWAStraining5_GPD.csv", index = None, header= 1)
traininggeno6.to_csv("GWAStraining6_GPD.csv", index = None, header= 1)






