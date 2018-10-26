#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB

import numpy as np
import pandas as pd
import random
import string


# Function to generate ID using characters and digits
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
np.savetxt("./toy_dataset/phenotype_toy_dataset.txt",  np.r_[[phenotype_header], phenotype_dataset], fmt='%s', delimiter=',')
