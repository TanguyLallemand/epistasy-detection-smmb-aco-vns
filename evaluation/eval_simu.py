#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2 BB
#          Jonathan Cruard M2 BB

import glob
import os
import numpy as np

###############################################################################
# This function permit to define some arguments and generate associated help
###############################################################################


def get_arguments():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", help="Give an input directory", type=str, action='store', required=True)
    parser.add_argument(
        "-o", "--output", help="Give an output directory", type=str, action='store', required=True)
    parser.add_argument(
        "-n", "--nruns", help="Number of method executions to be performed on each file in the dataset", type=int)
    args = parser.parse_args()
    return args

###############################################################################
# This function permit to get all txt file containned in input directory
###############################################################################


def get_input_files(input_directory):
    # Search for file ending with txt extension in a given directory
    input_files = glob.glob(input_directory + '*.txt')
    return input_files

###############################################################################
# This function permit to determine for every results if it is a FP, TP or FN
###############################################################################


def parsing_result_file(result_file, pattern):
    result = []
    false_negative = 0
    true_positive = 0
    false_positive = 0
    for line in result_file:
        if line.strip() == pattern:
            result.append('TP')
            true_positive += 1
        elif line.strip() == '' or len(line.strip()) != len(pattern):
            result.append('FN')
            false_negative += 1
        elif line.strip() != pattern:
            result.append('FP')
            false_positive += 1
    return [result, true_positive, false_negative, false_positive]

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
# This function permit to create an output file and writing results
###############################################################################


def creation_of_output_file(output_directory, results_file_parsed, name_of_file):
    base = os.path.basename(name_of_file)
    name = os.path.splitext(base)[0] + '_results.txt'
    if os.path.isdir(output_directory):
        name = output_directory + '/' + name
        with open(name, 'w') as output_file:
            for item in results_file_parsed:
                output_file.write(item + '\n')

###############################################################################
# This function permit to calcul recall
###############################################################################


def calc_recall(true_positive, false_negative):
    recall = true_positive / (true_positive + false_negative)
    return recall

###############################################################################
# This function permit to calcul precision
###############################################################################


def calc_precision(true_positive, false_positive):
    precision = true_positive / (true_positive + false_positive)
    return precision

###############################################################################
# This function permit to create measure file and calcul f_measure
###############################################################################


def creation_of_measure_file(recall, precision):
    with open('f_measures.txt', 'w') as measure_file:
        f_measure = 2 / (1 + recall + 1 + precision)
        measure_file.write(str(f_measure))
        return f_measure

###############################################################################
# This function permit to powers file and calcul powers
###############################################################################


def creation_of_powers_file(f_measure, number_of_execution, number_true_positive):
    with open('powers.txt', 'w') as powers_file:
        powers_file.write('Power\n')
        powers_file.write(
            str(number_true_positive / number_of_execution) + '\n')
        powers_file.write('f_measures\n')
        for item in f_measure:
            powers_file.write(str(item) + '\n')

###############################################################################
# Main function
###############################################################################


def main():
    # Initialization of some variables
    # Put pattern in raw code, TODO: get it
    pattern = 'ml'
    f_measure = []
    # Get argument parser
    args = get_arguments()
    # Get argument passed to script
    input_directory = args.input
    output_directory = args.output
    number_of_execution = args.nruns

    # Get list of file in input directory
    input_files = get_input_files(input_directory)
    # Get number of files in input directory
    n_files = len(input_files)

    # For every files
    for file in input_files:
        with open(file, 'r') as result_file:
            # Determine results category: TP, FN or FP
            results_file_parsed = parsing_result_file(result_file, pattern)
            # Check if output directory exist, if not creation of it
            check_output_directory(output_directory)
            # Save parsed result in a file
            creation_of_output_file(
                output_directory, results_file_parsed[0], file)

            # Save number of TP, FN and TN
            number_true_positive = results_file_parsed[1]
            number_false_negative = results_file_parsed[2]
            number_false_positive = results_file_parsed[3]
            # Resolve problems for calculating f-measure, if a variable is equal to zero add 1 and remove 1 to the other one
            if number_false_negative == 0:
                number_false_positive -= 1
                number_false_negative += 1
            if number_false_positive == 0:
                number_false_negative -= 1
                number_false_positive += 1

            # Calcul recall
            recall = calc_recall(number_true_positive, number_false_negative)
            # Calcul precision
            precision = calc_precision(
                number_true_positive, number_false_positive)
            # Stock in a list all f_measures calculated
            f_measure.append(creation_of_measure_file(recall, precision))
            # Creation of power file containning all f_measure, and calculated power
            creation_of_powers_file(
                f_measure, number_of_execution, number_true_positive)


if __name__ == "__main__":
    main()
