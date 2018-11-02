#!/usr/bin/env python3
# Does not need any shebang because used in a virtual environment

# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB

import numpy as np
import pandas as pd
import random
import string


###############################################################################
# Function to define arguments for script
###############################################################################

# required=True
def get_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--output", help="Give an output directory", type=str, action='store')
    parser.add_argument(
        "-p", "--prefix", help="Common prefix for files", type=str, action='store')
    parser.add_argument(
        "-f", "--file", help="Number of file to generate", type=int)
    parser.add_argument(
        "-v", "--variable", help="Number of variable to generate", type=int)
    parser.add_argument(
        "-pa", "--case", help="Number of case to generate", type=int)
    parser.add_argument(
        "-c", "--control", help="Number of control to generate", type=int)
    args = parser.parse_args()
    return args

###############################################################################
# Function to generate ID using characters and digits
###############################################################################


def generate_id(common_prefix, i):
    id = ''.join(common_prefix + str(i))
    return id

###############################################################################
# This function permit to create an output directory if it does not exist
###############################################################################


def check_output_directory(output_directory):
    import re
    import os
    # Import errno to handle with errors during directory creation
    from errno import EEXIST
    # Get current path and add sub directory name
    # Get current directory path
    my_path = os.getcwd() + '/' + output_directory
    # Try to create a new directory
    try:
        # Make a directory following path given
        os.mkdir(my_path)
    except OSError as exc:
        if exc.errno == EEXIST:
            pass
        else:
            raise

###############################################################################
# Generation of genotype dataset
###############################################################################


def generate_genotype_dataset(output_directory, id, number_of_variable, number_of_patient):

    # intitialisation array of ID
    genotype_id = []
    check_output_directory(output_directory)
    path = './' + output_directory + '/genotype_toy_dataset' + id + '.txt'
    #Generate id of variables
    for j in range(0, number_of_variable):
        genotype_id.append(''.join('SNP' + str(j)))

    # Generate dataset
    genotype_dataset = np.random.randint(
        3, size=(number_of_patient, number_of_variable))

    # Concatenate both arrays and print it as a csv file called genotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[[genotype_id], genotype_dataset], fmt='%s', delimiter=',')


###############################################################################
# Generation of phenotype dataset
###############################################################################


def generate_phenotype_dataset():

    # Generate header
    phenotype_header = ["Class"]
    # Generate dataset
    phenotype_dataset = (np.random.randint(2, size=(4000, 1)))
    # Concatenate both arrays and print it as a csv file called phenotype_toy_dataset.txt
    np.savetxt("./toy_dataset/phenotype_toy_dataset.txt",
               np.r_[[phenotype_header], phenotype_dataset], fmt='%s', delimiter=',')


###############################################################################
# Main function
###############################################################################


def main():
    # Get variables from arguments
    args = get_arguments()
    # Get argument passed to script
    output_directory = args.output
    common_prefix = args.prefix
    number_of_file = args.file
    number_of_variable = args.variable
    number_of_case = args.case
    number_of_control = args.control
    for i in range(0, number_of_file):
        id = generate_id(common_prefix, i)
        generate_genotype_dataset(
            output_directory, id, number_of_variable, number_of_case)


if __name__ == "__main__":
    main()

# Votre logiciel de simulation simu_naive.xxx prendra en entrée :
# un nom de répertoire de sortie
# le nombre de fichiers à générer
# le préfixe commun aux identifiants des fichiers
# le nombre de variables à simuler
# le nombre de cas à simuler
# le nombre de contrôles à simuler
# Vous fournirez un script Python launch_simu_naiv
