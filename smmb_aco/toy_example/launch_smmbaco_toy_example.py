#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
import os
import subprocess

# This programm require two parameters:
# First: path of genotype dataset
# Second: path of associated phenotype dataset
os.system('./smmb_aco/smmb_aco.exe ./smmb_aco/toy_example/toy_dataset/genotypes_toy_dataset.txt ./smmb_aco/toy_example/toy_dataset/phenotypes_toy_dataset.txt ./smmb_aco/toy_example/parameters.txt')
