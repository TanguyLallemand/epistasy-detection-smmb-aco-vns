#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
import os

# This programm require three parameters:
# First: path of genotype data set
# Second: path of associated phenotype data set
# Third: path to parameters file
os.system('./vns/vns.exe ./vns/toy_example/toy_dataset/genotypes_toy_dataset.txt ./vns/toy_example/toy_dataset/phenotypes_toy_dataset.txt ./vns/toy_example/parameters.txt')
