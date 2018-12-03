#!/usr/bin/env python3
import argparse
import numpy
import string
import random
from math import exp


###################################################################################################
###################################################################################################


##########
# Parser #
##########


argParser=argparse.ArgumentParser()
# Output directory
argParser.add_argument("-o", "--output",
					required=True,
					help="Output results directory",
					dest="output")
# Number of generated simulation
argParser.add_argument("-n", "--number",
					required=True,
					type=int,
					help="Number of generation to simulate",
					dest="number_sim")
# Common prefix files
argParser.add_argument("-p", "--prefix",
					required=True,
					help="Common prefix for each result file",
					dest="prefix")
# Number of simulated variables
argParser.add_argument("-v", "--variables",
					required=True,
					type=int,
					help="Total number of SNP",
					dest="variables")
# Number of simulated cases
argParser.add_argument("-c", "--cases",
					required=True,
					type=int,
					help="Number of simulated cases",
					dest="cases")
# Number of simulated controls
argParser.add_argument("-t", "--controls",
					required=True,
					type=int,
					help="Number of simulated controls",
					dest="controls")
# Size of epistasic pattern
argParser.add_argument("-s", "--size",
					required=True,
					type=int,
					help="Size of the epistasic pattern",
					dest="pattern_size")
# Args variable for every dest of the parse
args=argParser.parse_args()


###################################################################################################
###################################################################################################


#############################
# Random matrice generation #
#############################


# Loop for the number of simulation
for nbr_sim in range(0,args.number_sim):
	nbr_sim+=1

	# Random SNP id generation
	variables=args.variables-args.pattern_size
	SNP_list=[]
	for i in range(0,variables):
		snp=""
		letter="N"
		digit=string.digits
		digit_1=random.choice(digit)
		digit_2=random.choice(digit)
		snp=letter+digit_1+digit_2
		SNP_list.append(snp)
	random_SNP=numpy.asarray(SNP_list)

	# Random genotypes generation
	total=args.cases+args.controls
	random_nbr=numpy.random.randint(3, size=(total,variables))

	# Random genotypes matrix formation
	matrix_random_genotypes=numpy.vstack([random_SNP,random_nbr])

	# Linear regression for phenotypes generation (for the two last genotypes columns)
	mean=0
	compteur=0
	while mean!=0.5:
		compteur+=1
		pheno_list=[]
		geno_dict={}
		for i in range(0,total):
			beta_val_list=[-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-0.9,-1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
			x_val_list=[0,1,2]
			x_list=[]
			beta_list=[]
			for j in range(args.pattern_size):
				k=random.choice(x_val_list)
				x_list.append(k)
			for j in range(args.pattern_size):
				k=random.choice(beta_val_list)
				beta_list.append(k)
			y_list=[]
			for j in range(0,args.pattern_size):
				y=(beta_list[j])*(x_list[j])
				y_list.append(y)
			y=sum(y_list)
			p=1/(1+exp(-y))
			if p>=0.5:
				pheno_list.append(1)
			elif p<0.5:
				pheno_list.append(0)
			geno_dict[i]=x_list
		mean=(sum(pheno_list)/total)
	
	# Recovery of the SNP genotypes in a numpy matrice
	regression_genotypes=numpy.array(list(geno_dict.values()))
	
	# Regression SNP id generation
	SNP_list=[]
	for i in range(0,args.pattern_size):
		snp=""
		letter="M"
		digit=string.digits
		digit_1=random.choice(digit)
		digit_2=random.choice(digit)
		snp=letter+digit_1+digit_2
		SNP_list.append(snp)
	regression_SNP=numpy.asarray(SNP_list)
	
	# Regression genotypes matrix formation
	matrix_regression_genotypes=numpy.vstack([regression_SNP,regression_genotypes])	

	# Add of the matrix_regression_genotypes to matrix_random_genotypes
	matrix_genotypes=numpy.column_stack([matrix_random_genotypes,matrix_regression_genotypes])

	# Phenotypes matrix formation
	pheno_list.insert(0,"#### 0 : sain ; 1 : malade ####")
	phenotypes=numpy.asarray(pheno_list)


###################################################################################################
###################################################################################################


##############################
# SMMB-ACO files generations #
##############################

	# Numpy SMMB-ACO file generation
	numpy.savetxt(args.output+'/{0}genotypes_sim_{1}.txt'.format(args.prefix,nbr_sim),
				matrix_genotypes,
				fmt='%s',
				delimiter=',')
	numpy.savetxt(args.output+'/{0}phenotypes_sim_{1}.txt'.format(args.prefix,nbr_sim),
				phenotypes,
				fmt='%s')
