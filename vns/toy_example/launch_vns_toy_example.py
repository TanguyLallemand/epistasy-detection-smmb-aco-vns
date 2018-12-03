#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
import os
import subprocess

# This programm require two parameters:
# First: path of genotype dataset
# Second: path of associated phenotype dataset
os.system('./vns/vns.exe ./vns/toy_example/toy_dataset/genotypes_toy_dataset.txt ./vns/toy_example/toy_dataset/phenotypes_toy_dataset.txt ./vns/toy_example/parameters.txt')
