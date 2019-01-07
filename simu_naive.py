#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB

import numpy as np
import pandas as pd
import string
import itertools
import re
import os
# Import errno to handle with errors during directory creation
from errno import EEXIST
from math import exp
from random import *


################################################################################
# Function to define arguments for script                                      #
# Generate associated help                                                     #
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
# This function allows to create an output directory if it does not exist      #
################################################################################


def check_output_directory(output_directory):
    # Get current path and add sub directory name
    # Get current directory path
    my_path = os.getcwd() + '/' + output_directory
    # Try to create a new directory
    try:
        # Make a directory following path given
        os.mkdir(my_path)
    #  Allow to handle with errors during directory creation
    except OSError as exc:
        if exc.errno == EEXIST:
            pass
        else:
            raise

################################################################################
# determine_treshold allows to return a list of th 10 percent threshold        #
################################################################################


def determine_treshold(all_combinations):
    threshold = ()
    percentage = int(len(all_combinations) * 0.1)
    # Fix percentage to 100 if it's inferior to it
    if percentage < 1:
        percentage = 1
    for x in range(-percentage, percentage):
        #  Compute threshold
        threshold = threshold + (int(len(all_combinations) / 2 + x),)
    return threshold

################################################################################
# randrange_float allows to generate random int in a particular range and      #
# a given pitch                                                                #
################################################################################


def randrange_float(start, stop, pitch):
    return randint(0, int((stop - start) / pitch)) * pitch + start

################################################################################
# fit_logit allows to generate an array of array of beta coeffecient           #
################################################################################


def fit_logit(pattern_size, all_combinations, threshold):
    max_iterations = 1000
    array_psi_global = []
    array_random_global = []
    # Check for reliability of generated combination
    if len(all_combinations) != pow(3, pattern_size):
        print("Error in combinations")
    for logit_iter in range(0, max_iterations):
        count_healthy = 0
        array_of_psi = []
        array_random = []
        #  Generate an array of pattern size of random between 0 and 1
        for i in range(0, pattern_size):
            array_random.append(randrange_float(0, 1, 0.1))
        # For every possible cases
        for i in range(0, pow(3, pattern_size)):
            # COmpute precision of each cases associated with random
            precision = compute_logit(array_random, all_combinations[i])
            #  If precision < 0.5 associate with healthy case (control)
            if precision < 0.5:
                count_healthy += 1
        if count_healthy in threshold:
            # If number of healty is in accepted values save it
            array_random_global.append(array_random)
    return array_random_global

################################################################################
# compute_logit allows to compute a logistique regression using possible values#
# and possible genotypes                                                       #
################################################################################


def compute_logit(list_random, combination):
    Y = -1
    for i in range(0, len(list_random)):
        Y = Y + list_random[i] * combination[i]
        if i != 0:
            multiplicate_Bs = multiplicate_Bs * list_random[i]
            multiplicate_Xs = multiplicate_Xs * combination[i]
        else:
            multiplicate_Bs = list_random[0]
            multiplicate_Xs = combination[0]
    # Compute logistic regression
    Y = Y + (multiplicate_Bs * multiplicate_Xs)
    # Compute associated precision
    precision = (1 / (1 + exp(-Y)))
    return precision

################################################################################
# generate_SNP_name allows to generate an array of non causal and causal SNPs  #
# names                                                                        #
################################################################################


def generate_SNP_name(number_of_variables, pattern_size):
    id = []
    # Generate id of variables
    for j in range(0, number_of_variables - pattern_size):
        # Join "SNP" string with iterator
        id.append(''.join('N' + str(j)))
    for j in range(0, pattern_size):
        # Join "SNP" string with iterator
        id.append(''.join('M' + str(j)))
    return id

################################################################################
# save_phenotype_dataset use phenotype dataset, contruct a filename and save   #
# it                                                                           #
################################################################################


def save_phenotype_dataset(i, output_directory, common_prefix, phenotype_dataset):
    # Transtyp phenotype dataset
    final_matrix_phenotype_dataset = np.asarray(phenotype_dataset)
    # Generate header
    header = ["Class"]
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/' + \
        str(i) + '_' + common_prefix + '_phenotype_toy_dataset' + '.txt'
    # Print it as a csv file called phenotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[header, final_matrix_phenotype_dataset], fmt='%s', delimiter=',')

################################################################################
# save_genotype_dataset, use header and matrix, merge them, contruct a file    #
# name and save genotype dataset                                               #
################################################################################


def save_genotype_dataset(i, output_directory, common_prefix, matrix_genotype_ID, matrix_ready_save):
    # Generate header
    # Check if output directory exist
    check_output_directory(output_directory)
    # Generate path to save txt file
    path = './' + output_directory + '/' + \
        str(i) + '_' + common_prefix + '_genotype_toy_dataset' + '.txt'
    # Print it as a csv file called phenotype_toy_dataset.txt
    np.savetxt(path,
               np.r_[[matrix_genotype_ID], matrix_ready_save], fmt='%s', delimiter=',')

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

    # List of different genotype
    list_value_genotype = [[0, 1, 2]]
    #  Generate all possible pattern
    for i in range(1, size_pattern):
        list_value_genotype.append([0, 1, 2])
    all_combinations = list(itertools.product(*list_value_genotype))
    print("######################## Generate possible combinations ########################")
    threshold = determine_treshold(all_combinations)
    print("############################### Compute threshold ##############################")
    logit = fit_logit(size_pattern, all_combinations, threshold)
    print("######################### Compute logistic regression ##########################")
    # This loop will determine the number of file for a run(with the same logit model)
    for i in range(1, number_of_file + 1):
        print("############################## Generate data set ###############################")
        # SNPs IDs matrix initialization
        matrix_genotype_non_causal_ID = []
        # Pattern matrix initialization
        pattern_genotype_case = []
        pattern_genotype_control = []
        # Phenotypes matrix initialization
        matrix_phenotype_case = []
        matrix_phenotype_control = []
        # Generate ID of SNPs
        matrix_genotype_ID = generate_SNP_name(
            number_of_variable, size_pattern)

        # Matrix random creation
        matrix_case_geno = np.random.randint(low=0, high=3, size=(
            number_of_case, number_of_variable - size_pattern), dtype="int")
        matrix_control_geno = np.random.randint(low=0, high=3, size=(
            number_of_control, number_of_variable - size_pattern), dtype="int")

        # Select first list of Betas
        list_random = logit[0]
        # Generation of geno for pattern and phenotype
        # while there is less genotype generated than the number of case or the number of control
        while ((len(pattern_genotype_case) < number_of_case) or (len(pattern_genotype_control) < number_of_control)):
            temp_matrix = np.random.randint(
                low=0, high=3, size=(1, size_pattern))
            pr_value = compute_logit(list_random, temp_matrix[0])
            # If there pr is upper 0.5, the generated individual is sick
            if pr_value > 0.5:
                if len(pattern_genotype_case) >= number_of_case:
                    continue
                pattern_genotype_case.append(temp_matrix[0])
                matrix_phenotype_case.append(1)
            else:
                if len(pattern_genotype_control) >= number_of_control:
                    continue
                pattern_genotype_control.append(temp_matrix[0])
                matrix_phenotype_control.append(0)

        # hstack will concatenate by column
        matrix_final_geno_case = np.hstack(
            (matrix_case_geno, pattern_genotype_case))
        matrix_final_geno_control = np.hstack(
            (matrix_control_geno, pattern_genotype_control))
        # vstack will concatenate by rows
        matrix_ready_save = np.vstack(
            (matrix_final_geno_case, matrix_final_geno_control))
        # concatenate pheno matrix
        matrix_final_pheno = np.hstack(
            (matrix_phenotype_case, matrix_phenotype_control))
        save_phenotype_dataset(i,
                               output_directory, common_prefix, matrix_final_pheno)
        save_genotype_dataset(i, output_directory, common_prefix,
                              matrix_genotype_ID, matrix_ready_save)


################################################################################
# Execution of main function                                                   #
################################################################################

if __name__ == "__main__":
    main()
