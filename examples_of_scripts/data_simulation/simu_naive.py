#! /usr/bin/python3
#-*-coding: utf-8-*-
# In order to launch the script in a conda environment, use this call: " Python ./file_to_launch"
# Louison Fresnais M2BB
# François Courtin M2BB

import sys
import re
import os


from parameter_parsing import *
from logistic_model import *
import pandas as pd
import numpy as np
args = vars(argument_parsing())
print(args)


#####################################################################
#                           INPUT                                     #
#####################################################################

###Number of file reading
try:
    number_file = int(args.get("file"))
    print("Number of file value: ",number_file)

except TypeError:
    number_file = int(input("enter the number of file you want to generate in one run"))

###Prefix value reading
prefix = str(args.get('prefix'))
print("Prefix value: ",prefix)
if (prefix == 'None'):

    prefix = str(input("enter a disctinctive prefix that will be applied to each result file for this run"))


###Number of variables reading
try:
    number_variable = int(args.get("variable"))
    print("Number of variable value: ",number_file)

except TypeError:
    number_variable = int(input("enter the number of variables to simulate"))

###Number of cases reading
try:
    number_case = int(args.get("case"))
    print("Number of case value: ",number_case)

except TypeError:
    number_case = int(input("enter the number of cases to simulate"))

###Number of controls reading
try:
    number_control = int(args.get("control"))
    print("Number of control value: ",number_control)

except TypeError:
    number_control = int(input("enter the number of controls to simulate"))

###Size of epistasia pattern
try:
    size_epistasia = int(args.get("pattern"))
    print("Size of the epistasia pattern: ",size_epistasia)
    if number_variable <= size_epistasia:
        number_variable = int(input("the pattern size must be greater than the number of variable, enter a new variable number  "))

except TypeError:
    size_epistasia = int(input("enter the size of the epistasia pattern"))
    if number_variable <= size_epistasia:
        number_variable = int(input("the pattern size must be greater than the number of variable, enter a new variable number "))

#####################################################################
#                           EXECUTION                                 #
#####################################################################
###Logistic Regression
#list of different value for a gene
list_comb = [[0,1,2]]

max_iter = 1000 #This value determine the number of iteration allowed to find a good logit model
#Generate all combinations according to list comb and size_epistasia
all_combinations = epistasis_combination(size_epistasia,list_comb)

#Determine the 10 percent threshold that will decide if there
#is a good  enough distribution of case/control with the logit model
determination_th = determination_th(all_combinations,0.1)

#Will provide a list of list of betas for the logit model
logit_list = fit_relevant_logit(size_epistasia,all_combinations,determination_th,max_iter)


#Folder creation with a bash command
try:
    folder = str(prefix+"_all_files")
    os.mkdir(folder)
except FileExistsError:
        folder = str(prefix+"_all_files"+str(np.random.randint(low=1,high=10000,size=(1))))
        os.mkdir(folder)


#This loop will determine the number of file for a run(with the same logit model)
for i in range(1,number_file+1):
    #file names creation
    geno_file_name = str(prefix+"genotypes"+str(i)+".txt")
    pheno_file_name = str(prefix+"phenotypes"+str(i)+".txt")
    #SNPs IDs matrix initialization
    matrix_genotype_ID = []
    #pattern matrix initialization
    pattern_genotype_case = []
    pattern_genotype_control = []
    #phenotypes matrix initialization
    matrix_phenotype_case = []
    matrix_phenotype_control = []

    #Creation of SNPs name's
    for x in range(1,number_variable+1):
        x = str(x)
        ID = "N"+x
        matrix_genotype_ID.append(ID)

    #Convert the ID list into an array
    matrix_genotype_ID = np.asarray(matrix_genotype_ID)
    #Matrix random creation
    matrix_case_geno = np.random.randint(low=0,high=3, size=(number_case, number_variable-size_epistasia),dtype="int")
    matrix_control_geno = np.random.randint(low=0,high=3, size=(number_control, number_variable-size_epistasia))

    #Select first list of Betas
    list_random = logit_list[1]

    #Generation of geno for pattern and phenotype
    #while there is less genotype generated than the number of case or the number of control
    while ((len(pattern_genotype_case) < number_case) or (len(pattern_genotype_control) < number_control)):
        temp_matrix = np.random.randint(low=0,high=3,size=(1,size_epistasia))
        pr_value = compute_logit(list_random,temp_matrix[0])
        #If there pr is upper 0.5, the generated individual is sick
        if pr_value > 0.5:
            if len(pattern_genotype_case) >= number_case:
                continue;
            pattern_genotype_case.append(temp_matrix[0])
            matrix_phenotype_case.append(1)
        else:
            if len(pattern_genotype_control) >= number_control:
                continue;
            pattern_genotype_control.append(temp_matrix[0])
            matrix_phenotype_control.append(0)

    #hstack will concatenate by column
    matrix_final_geno_case = np.hstack((matrix_case_geno,pattern_genotype_case))
    matrix_final_geno_control = np.hstack((matrix_control_geno,pattern_genotype_control))
    #vstack will concatenate by rows
    matrix_ready_save = np.vstack((matrix_final_geno_case,matrix_final_geno_control))

    #save the simulated genotype data file
    np.savetxt(folder+"/"+geno_file_name,np.r_[[matrix_genotype_ID], matrix_ready_save], fmt='%s', delimiter=',')
    header_pheno_file=["Class"]
    #concatenate pheno matrix
    matrix_final_pheno = np.hstack((matrix_phenotype_case,matrix_phenotype_control))
    #save the simulated phenotype data file
    np.savetxt(folder+"/"+pheno_file_name,np.r_[header_pheno_file,matrix_final_pheno], fmt='%s', delimiter=',')
