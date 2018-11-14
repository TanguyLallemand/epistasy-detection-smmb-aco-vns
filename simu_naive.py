#!/usr/bin/env python3

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


def get_arguments():
    import argparse
    # Define arument parser and stock it as object
    parser = argparse.ArgumentParser()
    # Add all arguments in parser
    parser.add_argument(
        "-o", "--output", help="Give an output directory", type=str, action='store', required=True)
    parser.add_argument(
        "-p", "--prefix", help="Common prefix for files", type=str, action='store', required=True)
    parser.add_argument(
        "-f", "--file", help="Number of file to generate", type=int, required=True)
    parser.add_argument(
        "-v", "--variable", help="Number of variable to generate", type=int, required=True)
    parser.add_argument(
        "-pa", "--case", help="Number of case to generate", type=int, required=True)
    parser.add_argument(
        "-c", "--control", help="Number of control to generate", type=int, required=True)
    # Parse arguments
    args = parser.parse_args()
    # Return parsed args
    return args

###############################################################################
# Function to generate ID using characters and digits
###############################################################################


def generate_id(common_prefix, i):
    # Generate an ID by joining common prefix and an iterator
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

    # Intitialisation array of ID
    genotype_id = []
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/genotype_toy_dataset' + id + '.txt'
    # Generate id of variables
    for j in range(0, number_of_variable):
        # Join "SNP" string with iterator
        genotype_id.append(''.join('SNP' + str(j)))

    # Generate dataset
    genotype_dataset = np.random.randint(
        3, size=(number_of_patient, number_of_variable))

    # Concatenate both arrays and print it as a csv file called genotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[[genotype_id], genotype_dataset], fmt='%s', delimiter=',')


###############################################################################
# Generation of positive SNPs
###############################################################################
def generate_causals_SNPs(pattern_size, ):



###############################################################################
# Generation of phenotype dataset
###############################################################################


def generate_phenotype_dataset(output_directory, number_of_case, number_of_control, common_prefix):
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/phenotype_toy_dataset' + common_prefix + '.txt'
    # Generate header
    phenotype_header = ["Class"]

    # Generate dataset
    # Generate right number of case as 1
    phenotype_case = np.ones((number_of_case, 1), dtype=int)
    # Generate right number of control as 0
    phenotype_control = np.zeros((number_of_control, 1), dtype=int)
    # Merge those two numpy array
    phenotype_dataset = np.append(phenotype_control, phenotype_case, axis=0)

    # Concatenate both arrays and print it as a csv file called phenotype_toy_dataset.txt
    np.savetxt(path,
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
    # Calculate number of patient
    number_of_patients = number_of_case + number_of_control
    # A loop to generate number of genotype file asked
    for i in range(0, number_of_file):
        # Generate an ID for current file in construction
        id = generate_id(common_prefix, i)
        # Generate datas and save it in a CSV file
        generate_genotype_dataset(
            output_directory, id, number_of_variable, number_of_patients)
    # Generate one phenotype dataset linked to generated dataset
    generate_phenotype_dataset(
        output_directory, number_of_case, number_of_control, common_prefix)


# Execution of main function
if __name__ == "__main__":
    main()
