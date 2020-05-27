#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gwf workflow found in source activate gwfenviron.

#Change only grouping1,grouping2,grouping3 to change round


This is a workflow used for genomic prediction
by first performing EMMAX+EMMA200 GWAS on a n-20 training 
population of clover, leaving 20 individuals out for testing,
then using the top 200 significant SNPs (MAF=0.05 filter)
to make a GRM used in GBLUP to predict testing population or 
in EpiG+GBLUP to predict GEBVs of testing population individuals.

That workflow consists of four parts.
Part1: Creating a training population and a testing population.
This part outputs 1) A GWAS-ready genotype with only training individuals

Part2: Running EMMAX+EMMA200 GWAS.
A .pval file is generated with p-values of all SNPs.

Part3: Making a GRM based on top 200 significant SNPs.
Outputs 1) a list of top 200 significant SNPs,
2) A GRM including all n individuals, 
based on top 200 significant SNPs only.

Part4: Genomic prediction using GBLUP
Genomic prediction is performed using the R-package
BGLR. 

Part5: All predictions are combined into one file to create one correlation coefficient
"""

from gwf import Workflow

gwf=Workflow()




def script_caller(groups, geno, outputs, script_name):
    inputs = [groups,geno]
    outputs = outputs
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '08:00:00'
    }

    spec = '''
    python {script_name} {groups} {geno}
    '''.format(script_name=script_name, groups=groups,geno=geno)

    return inputs, outputs, options, spec




def GRM_script_caller(pval, geno,top200,GRM):
    inputs = [pval,geno]
    outputs = [top200,GRM]
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '08:00:00'
    }

    spec = '''
    python Make_GRM_based_on_top200_SNP.py {pval} {geno} {top200} {GRM}
    '''.format(pval=pval,geno=geno,top200=top200,GRM=GRM)
    
    
    return inputs, outputs, options, spec



def GWAS_caller(phenofile, genofile,dicname,outputs):
    inputs = [phenofile,genofile]
    outputs = outputs
    options = {
        'cores': 9,
        'memory': '8g',
        'walltime': '08:00:00'
    }
    

    spec = '''
    source ~/miniconda2/etc/profile.d/conda.sh
    conda activate myproject
    mkdir {dicname}
    cd {dicname}
    python /home/cks/NChain/faststorage/GWASimplementation/atgwas/src/gwa.py -o "" -a emmax -m 6 -r {phenofile} -f {genofile} --data_format="diploid_int" 
    '''.format(dicname=dicname,phenofile=phenofile,genofile=genofile)
    
    return inputs, outputs, options, spec


def GP_script_caller(GRM,name,trainingpop,outputs):
    inputs= [GRM,name,trainingpop]
    outputs=outputs
    options = {
        'cores': 16,
        'memory': '100g',
        'walltime': '02:00:00'    
    }
    
    spec='''
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate Rprogram
    Rscript GBLUP_replicate_data_GPD.R {GRM} {name} {trainingpop} 
    '''.format(GRM=GRM,name=name,trainingpop=trainingpop)
    
    return inputs, outputs, options, spec


def Correlation_calc(results1,results2,results3,results4,results5,results6,outputs):
    inputs = [results1,results2,results3,results4,results5,results6]
    outputs = outputs
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '02:00:00'    
    }

    spec='''
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate Rprogram
    Rscript Correlation.R {results1} {results2} {results3} {results4} {results5} {results6} 
    '''.format(results1=results1,results2=results2,results3=results3,results4=results4,results5=results5,results6=results6)
    
    return inputs, outputs, options, spec
    



gwf.target_from_template("Training_pop_generator",
                         script_caller(groups="../grouping1.txt",geno="RNAseq_LDfiltered_20190823.csv",
                                       outputs=["GWAStraining1_GPD.csv","GWAStraining2_GPD.csv","GWAStraining3_GPD.csv","GWAStraining4_GPD.csv","GWAStraining5_GPD.csv","GWAStraining6_GPD.csv"],
                                       script_name="GWAS_creating_trainingpop_20191002.py"))

gwf.target_from_template("GWAS_1",
                         GWAS_caller(dicname="round1_tstpop1",phenofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/gpd_dryweight_mean_cor_20190923.csv",
                         genofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/GWAStraining1_GPD.csv",
                         outputs=["/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop1/results/_pid1_GPD_corrected_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_2",
                         GWAS_caller(dicname="round1_tstpop2",phenofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/gpd_dryweight_mean_cor_20190923.csv",
                         genofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/GWAStraining2_GPD.csv",
                         outputs=["/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop2/results/_pid1_GPD_corrected_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_3",
                         GWAS_caller(dicname="round1_tstpop3",phenofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/gpd_dryweight_mean_cor_20190923.csv",
                         genofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/GWAStraining3_GPD.csv",
                         outputs=["/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop3/results/_pid1_GPD_corrected_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_4",
                         GWAS_caller(dicname="round1_tstpop4",phenofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/gpd_dryweight_mean_cor_20190923.csv",
                         genofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/GWAStraining4_GPD.csv",
                         outputs=["/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop4/results/_pid1_GPD_corrected_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_5",
                         GWAS_caller(dicname="round1_tstpop5",phenofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/gpd_dryweight_mean_cor_20190923.csv",
                         genofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/GWAStraining5_GPD.csv",
                         outputs=["/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop5/results/_pid1_GPD_corrected_emmax_none_t75.pvals"]))

gwf.target_from_template("GWAS_6",
                         GWAS_caller(dicname="round1_tstpop6",phenofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/gpd_dryweight_mean_cor_20190923.csv",
                         genofile="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/GWAStraining6_GPD.csv",
                         outputs=["/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop6/results/_pid1_GPD_corrected_emmax_none_t75.pvals"]))

gwf.target_from_template("GRM1",
                         GRM_script_caller(pval="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop1/results/_pid1_GPD_corrected_emmax_none_t75.pvals",geno="RNAseq_LDfiltered_20190823.csv",
                                           top200="Round1_tstpop1_top200snps.txt",GRM='GRMbasedontop200SNPs_Round1_tstpop1.csv'))

gwf.target_from_template("GRM2",
                         GRM_script_caller(pval="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop2/results/_pid1_GPD_corrected_emmax_none_t75.pvals",geno="RNAseq_LDfiltered_20190823.csv",
                                           top200="Round1_tstpop2_top200snps.txt",GRM='GRMbasedontop200SNPs_Round1_tstpop2.csv'))

gwf.target_from_template("GRM3",
                         GRM_script_caller(pval="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop3/results/_pid1_GPD_corrected_emmax_none_t75.pvals",geno="RNAseq_LDfiltered_20190823.csv",
                                           top200="Round1_tstpop3_top200snps.txt",GRM='GRMbasedontop200SNPs_Round1_tstpop3.csv'))

gwf.target_from_template("GRM4",
                         GRM_script_caller(pval="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop4/results/_pid1_GPD_corrected_emmax_none_t75.pvals",geno="RNAseq_LDfiltered_20190823.csv",
                                           top200="Round1_tstpop4_top200snps.txt",GRM='GRMbasedontop200SNPs_Round1_tstpop4.csv'))

gwf.target_from_template("GRM5",
                         GRM_script_caller(pval="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop5/results/_pid1_GPD_corrected_emmax_none_t75.pvals",geno="RNAseq_LDfiltered_20190823.csv",
                                           top200="Round1_tstpop5_top200snps.txt",GRM='GRMbasedontop200SNPs_Round1_tstpop5.csv'))

gwf.target_from_template("GRM6",
                         GRM_script_caller(pval="/faststorage/project/NChain/WHITE_CLOVER/RNASEQ/HyLiTE/LD_filtering_RNASeq/GP_20190919_GPD_onReplicateData/GWASbased_GP_20191002/round1_tstpop6/results/_pid1_GPD_corrected_emmax_none_t75.pvals",geno="RNAseq_LDfiltered_20190823.csv",
                                           top200="Round1_tstpop6_top200snps.txt",GRM='GRMbasedontop200SNPs_Round1_tstpop6.csv'))

gwf.target_from_template("GP_round1tstpop1",
                         GP_script_caller(GRM="GRMbasedontop200SNPs_Round1_tstpop1.csv",
                                          name="round1_tstpop1",
                                          trainingpop="GWAStraining1_GPD.csv",
                                          outputs=["Predictions_GBLUP_GPDround1_tstpop1.txt"]))

gwf.target_from_template("GP_round1tstpop2",
                         GP_script_caller(GRM="GRMbasedontop200SNPs_Round1_tstpop2.csv",
                                          name="round1_tstpop2",
                                          trainingpop="GWAStraining2_GPD.csv",
                                          outputs=["Predictions_GBLUP_GPDround1_tstpop2.txt"]))

gwf.target_from_template("GP_round1tstpop3",
                         GP_script_caller(GRM="GRMbasedontop200SNPs_Round1_tstpop3.csv",
                                          name="round1_tstpop3",
                                          trainingpop="GWAStraining3_GPD.csv",
                                          outputs=["Predictions_GBLUP_GPDround1_tstpop3.txt"]))

gwf.target_from_template("GP_round1tstpop4",
                         GP_script_caller(GRM="GRMbasedontop200SNPs_Round1_tstpop4.csv",
                                          name="round1_tstpop4",
                                          trainingpop="GWAStraining4_GPD.csv",
                                          outputs=["Predictions_GBLUP_GPDround1_tstpop4.txt"]))

gwf.target_from_template("GP_round1tstpop5",
                         GP_script_caller(GRM="GRMbasedontop200SNPs_Round1_tstpop5.csv",
                                          name="round1_tstpop5",
                                          trainingpop="GWAStraining5_GPD.csv",
                                          outputs=["Predictions_GBLUP_GPDround1_tstpop5.txt"]))

gwf.target_from_template("GP_round1tstpop6",
                         GP_script_caller(GRM="GRMbasedontop200SNPs_Round1_tstpop6.csv",
                                          name="round1_tstpop6",
                                          trainingpop="GWAStraining6_GPD.csv",
                                          outputs=["Predictions_GBLUP_GPDround1_tstpop6.txt"]))

gwf.target_from_template("Cor_round3",
                         Correlation_calc(results1="Predictions_GBLUP_GPDround1_tstpop1.txt",results2="Predictions_GBLUP_GPDround1_tstpop2.txt",results3="Predictions_GBLUP_GPDround1_tstpop3.txt",results4="Predictions_GBLUP_GPDround1_tstpop4.txt",results5="Predictions_GBLUP_GPDround1_tstpop5.txt",results6="Predictions_GBLUP_GPDround1_tstpop6.txt",
                                          outputs=["Predictions_GBLUP_GPD1.txt","Correlation1.txt"]))
















