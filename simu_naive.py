# Votre logiciel de simulation simu_naive.xxx prendra en entrée :
# un nom de répertoire de sortie
# le nombre de fichiers à générer
# le préfixe commun aux identifiants des fichiers
# le nombre de variables à simuler
# le nombre de cas à simuler
# le nombre de contrôles à simuler
# Vous fournirez un script Python launch_simu_naiv

#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB

import numpy as np
import pandas as pd
import random
import string
import argparse

def get_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--output", help="Give an output directory", type=str, action='store', required=True)
    parser.add_argument(
        "-p", "--prefix", help="Common prefix for files", type=str, action='store', required=True)
    parser.add_argument(
        "-n", "--file", help="Number of file to generate", type=int)
    parser.add_argument(
        "-v", "--variable", help="Number of variable to generate", type=int)
    parser.add_argument(
        "-p", "--patient", help="Number of patient to generate", type=int)
    parser.add_argument(
        "-c", "--control", help="Number of control to generate", type=int)
    args = parser.parse_args()
    return args

###############################################################################
# Function to generate ID using characters and digits
###############################################################################


def random_id(size):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(size))


###############################################################################
# Generation of genotype dataset
###############################################################################
# intitialisation array of ID
genotype_id = []
# Loop to get 28 ID
for i in range(28):
    genotype_id.append(random_id(2))
# Generate dataset
genotype_dataset = np.random.randint(3, size=(4000, 28))
# Concatenate both arrays and print it as a csv file called genotype_toy_dataset.txt
np.savetxt("./toy_dataset/genotype_toy_dataset.txt",
           np.r_[[genotype_id], genotype_dataset], fmt='%s', delimiter=',')

###############################################################################
# Generation of phenotype dataset
###############################################################################
# Generate header
phenotype_header = ["Class"]
# Generate dataset
phenotype_dataset = (np.random.randint(2, size=(4000, 1)))
# Concatenate both arrays and print it as a csv file called phenotype_toy_dataset.txt
np.savetxt("./toy_dataset/phenotype_toy_dataset.txt",
           np.r_[[phenotype_header], phenotype_dataset], fmt='%s', delimiter=',')
