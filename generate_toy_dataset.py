#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2 BB
#          Jonathan Cruard M BB

import numpy as np
import pandas as pd
import random
import string


def random_id(size):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(size))


genotype_id = []
for i in range(28):
    genotype_id.append(random_id(2))
genotype_dataset = np.random.randint(3, size=(4000, 28))
np.savetxt("./toy_dataset/genotype_toy_dataset.txt",genotype_dataset, fmt='%i')
np.savetxt("./toy_dataset/genotype_toy_dataset.txt", genotype_id, "%s")
