#!/usr/bin/env python3

# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB

import numpy as np
import pandas as pd
import random
import string
from math import exp


###############################################################################
# Function to define arguments for script
###############################################################################


def get_arguments():
    import argparse
    # Define arument parser and stock it as object
    parser = argparse.ArgumentParser()
    # Add all arguments in parser
    parser.add_argument(
        "-o", "--output",
        help="Give an output directory",
        type=str,
        action='store',
        required=True)
    parser.add_argument(
        "-p", "--prefix",
        help="Common prefix for files",
        type=str,
        action='store',
        required=True)
    parser.add_argument(
        "-f", "--file",
        help="Number of file to generate",
        type=int,
        required=True)
    parser.add_argument(
        "-v", "--variable",
        help="Number of variable to generate",
        type=int,
        required=True)
    parser.add_argument(
        "-pa", "--case",
        help="Number of case to generate",
        type=int,
        required=True)
    parser.add_argument(
        "-c", "--control",
        help="Number of control to generate",
        type=int,
        required=True)
    parser.add_argument(
        "-s", "--size_pattern",
        help="Size of epistasis pattern to generate",
        type=int,
        required=True)
    # Parse arguments
    args = parser.parse_args()
    # Return parsed args
    return args

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


# Linear regression for phenotypes generation (for the two last genotypes columns)
def generate_genotype_with_linear_regression(number_of_patient):
    mean = 0
    while mean != 0.5:
        pheno_list = []
        genotype_linearly_generated_dataset = {}
        for i in range(0, number_of_patient):
            array_of_beta = [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -
                             0.8, -0.9, -1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            possible_values = [0, 1, 2]
            x_list = []
            list_of_beta = []
            y_list = []
            for j in range(args.pattern_size):
                k = random.choice(possible_values)
                x_list.append(k)
            for j in range(args.pattern_size):
                k = random.choice(beta_val_list)
                list_of_beta.append(k)
            for j in range(0, args.pattern_size):
                y = (list_of_beta[j]) * (x_list[j])
                y_list.append(y)
            y = sum(y_list)
            p = 1 / (1 + exp(-y))
            if p >= 0.5:
                pheno_list.append(1)
            elif p < 0.5:
                pheno_list.append(0)
            genotype_linearly_generated_dataset[i] = x_list
        mean = (sum(pheno_list) / number_of_patient)

    return [list(genotype_linearly_generated_dataset.values()), pheno_list]
###############################################################################
# Generation of genotype dataset
###############################################################################


def generate_genotype_dataset(output_directory, id, number_of_variable, number_of_patient, pattern_size):

    # Intitialisation array of ID
    random_genotype_id = []
    causal_genotype_id = []
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/genotype_toy_dataset' + id + '.txt'
    non_causal = number_of_variable - pattern_size
    # Generate id of variables
    for j in range(0, number_of_variable - non_causal):
        # Join "SNP" string with iterator
        random_genotype_id.append(''.join('SNP-N-' + str(j)))

    for j in range(0, pattern_size):
        # Join "SNP" string with iterator
        causal_genotype_id.append(''.join('SNP-C-' + str(j)))
    # Generate dataset
    random_genotype_dataset = np.random.randint(
        3, size=(number_of_patient, non_causal))
    # Merge SNP header with random generated matrix
    matrix_random_genotypes = numpy.vstack(
        [random_genotype_id, genotype_dataset])
    # Generate last SNP which are causaux
    results = generate_genotype_with_linear_regression(number_of_patient)
    linear_regression_genotype_dataset=np.array(results[0])
    # Merge SNP header with linear generated matrix
    matrix_linear_genotypes = numpy.vstack(
        [causal_genotype_id, matrix_random_genotypes])
    # Merge both matrix into one
    genotype_dataset=numpy.column_stack([matrix_linear_genotypes,matrix_random_genotypes])
    # Concatenate both arrays and print it as a csv file called genotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[genotype_dataset], fmt='%s', delimiter=',')


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
    id_non_causal[]
    id_causaux[]
    # Get variables from arguments
    args = get_arguments()
    # Get argument passed to script
    output_directory = args.output
    common_prefix = args.prefix
    number_of_file = args.file
    number_of_variable = args.variable
    number_of_case = args.case
    number_of_control = args.control
    size_pattern = args.size_pattern
    # Calculate number of patient
    number_of_patients = number_of_case + number_of_control
    # A loop to generate number of genotype file asked
    for i in range(0, number_of_file):
        # Generate datas and save it in a CSV file
        generate_genotype_dataset(
            output_directory, id, number_of_variable, number_of_patients, size_pattern)
    # Generate one phenotype dataset linked to generated dataset
    generate_phenotype_dataset(
        output_directory, number_of_case, number_of_control, common_prefix)


# Execution of main function
if __name__ == "__main__":
    main()
