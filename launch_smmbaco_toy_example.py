import os
import subprocess

# This programm require two parameters:
# First: path of genotype dataset
# Second: path of associated phenotype dataset
subprocess.Popen('./smmb_aco/smmb_aco.exe ./toy_dataset/genotype_toy_dataset_simu_naive0.txt ./toy_dataset/phenotype_toy_dataset_simu_naive.txt', shell=True)
