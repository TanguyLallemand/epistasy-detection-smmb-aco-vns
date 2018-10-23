#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2 BB
#          Jonathan Cruard M BB

import numpy as np
import pandas as pd
import random
import string


# Function to generate ID using characters and digits
def random_id(size):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(size))


# intitialisation array of ID
genotype_id = []
# Loop to get 28 ID
for i in range(28):
    genotype_id.append(random_id(2))
# Generate dataset
genotype_dataset = np.random.randint(3, size=(4000, 28))

# Concatenate both arrays and print it as a csv file
np.savetxt("./toy_dataset/genotype_toy_dataset.txt",
           np.r_[[genotype_id], genotype_dataset], fmt='%s', delimiter=',')
