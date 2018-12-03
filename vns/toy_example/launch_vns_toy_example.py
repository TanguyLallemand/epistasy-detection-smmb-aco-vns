#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
import os
import subprocess

# This programm require two parameters:
# First: path of genotype dataset
# Second: path of associated phenotype dataset
subprocess.Popen('./vns/vns.exe ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt', shell=True)
