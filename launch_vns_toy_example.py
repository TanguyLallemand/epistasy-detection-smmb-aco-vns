import os
import subprocess

# This programm require two parameters:
# First: path of genotype dataset
# Second: path of associated phenotype dataset
subprocess.Popen('./vns/vns.exe ./toy_dataset/genotypes_toy_dataset.txt ./toy_dataset/phenotypes_toy_dataset.txt', shell=True)
