#!/usr/bin/env python3

# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB

import numpy as np
import pandas as pd
import random
import string
from math import exp


################################################################################
# Function to define arguments for script                                      #
################################################################################


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

################################################################################
# This function permit to create an output directory if it does not exist      #
################################################################################


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

################################################################################
# Generation of causal genotype dataset                                        #
################################################################################


def generate_genotype_with_linear_regression(number_of_patient, pattern_size):
    # Initialisation of some variables
    # Variable to store mean
    mean = 0
    # List of possible value for beta
    array_of_beta = [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -
                     0.8, -0.9, -1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    # List of possible value in genotype matrix: 0: absence of SNP,
    # 1: SNP heterozygote, 2: SNP homozygote
    possible_values = [0, 1, 2]
    print("############################ Generating causal SNPs ############################")
    # Loop for generate right values
    while mean != 0.5:
        # Initialisation of temporary lists
        associated_phenotype = []
        genotype_linearly_generated_dataset = {}
        # For every values to generate
        for i in range(0, number_of_patient):
            # Initialisation of list to store temporary values
            x_list = []
            list_of_beta = []
            y_list = []
            # For every causal SNPs to generate
            for j in range(pattern_size):
                # Generate a value to fill matrix
                x_list.append(random.choice(possible_values))
            for j in range(pattern_size):
                # Choose a beta
                list_of_beta.append(random.choice(array_of_beta))
            for j in range(0, pattern_size):
                # Calculate y of every generated tuple matrix value and associated beta
                y_list.append((list_of_beta[j]) * (x_list[j]))
            # Sum all y calculated
            y = sum(y_list)
            # Calcluate precision to determine phenotype
            precision = 1 / (1 + exp(-y))
            # If precision is greater than 0.5 patient is sick
            if precision >= 0.5:
                associated_phenotype.append(1)
            # If precision is lesser than 0.5 patient is a control
            elif precision < 0.5:
                associated_phenotype.append(0)
            # Stor generated values
            genotype_linearly_generated_dataset[i] = x_list
        #Calculate mean to get 50% of sick
        mean = (sum(associated_phenotype) / number_of_patient)
        # Return genotype dataset and associated phenotype
    return [list(genotype_linearly_generated_dataset.values()), associated_phenotype]

################################################################################
# Generation of genotype dataset                                               #
################################################################################


def generate_genotype_dataset(output_directory, number_of_variable, number_of_patients, common_prefix, pattern_size, causal_genotype_SNPs):
    # Initialisation of list to store SNPs IDs
    random_genotype_id = []
    causal_genotype_id = []
    # Get causal_genotype_SNPs and transtype it in array
    linear_regression_genotype_dataset = np.array(causal_genotype_SNPs)
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/genotype_toy_dataset' + common_prefix + '.txt'
    non_causal = number_of_variable - pattern_size
    # Generate id of variables
    for j in range(0, non_causal):
        # Join "SNP" string with iterator
        random_genotype_id.append(''.join('SNP-N-' + str(j)))
    for j in range(0, pattern_size):
        # Join "SNP" string with iterator
        causal_genotype_id.append(''.join('SNP-C-' + str(j)))
    # Generate dataset
    random_genotype_dataset = np.random.randint(
        3, size=(number_of_patients, non_causal))
    # Merge SNP header with random generated matrix
    matrix_random_genotypes = np.vstack(
        [random_genotype_id, random_genotype_dataset])
    # Merge SNP header with linear generated matrix
    matrix_linear_genotypes = np.vstack(
        [causal_genotype_id, linear_regression_genotype_dataset])
    # Merge both matrix into one
    genotype_dataset = np.column_stack(
        [matrix_random_genotypes, matrix_linear_genotypes])
    # Concatenate both arrays and print it as a csv file called genotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[genotype_dataset], fmt='%s', delimiter=',')


################################################################################
# Generation of phenotype dataset                                              #
################################################################################


def generate_phenotype_dataset(output_directory, number_of_case, number_of_control, common_prefix, causal_phenotype_SNPs):
    # Generate header
    causal_phenotype_SNPs.insert(0, "Class ###### 0: healthy 1: sick ######")
    final_matrix_causal_phenotype_SNPs = np.asarray(causal_phenotype_SNPs)
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/phenotype_toy_dataset' + common_prefix + '.txt'
    # Print it as a csv file called phenotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[final_matrix_causal_phenotype_SNPs], fmt='%s', delimiter=',')


################################################################################
# Main function                                                                #
################################################################################


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
    size_pattern = args.size_pattern
    # Calculate number of patient
    number_of_patients = number_of_case + number_of_control
    # A loop to generate number of genotype file asked
    for i in range(0, number_of_file):
        # Generate a genotype matrix with causal SNPs
        results_of_regression = generate_genotype_with_linear_regression(
            number_of_patients, size_pattern)
        # Parse results
        # Get genotype dataset
        causal_genotype_SNPs = results_of_regression[0]
        # Get phenotype dataset
        causal_phenotype_SNPs = results_of_regression[1]
        # Generate datas and save it in a CSV file
        generate_genotype_dataset(
            output_directory, number_of_variable, number_of_patients, common_prefix, size_pattern, causal_genotype_SNPs)
        # Generate one phenotype dataset linked to generated dataset
        generate_phenotype_dataset(
            output_directory, number_of_case, number_of_control, common_prefix, causal_phenotype_SNPs)


################################################################################
# Execution of main function                                                   #
################################################################################
if __name__ == "__main__":
    main()
