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
    parser.add_argument(
        "-m", "--method", help="Method to test", type=str, action='store', required=True)
    args = parser.parse_args()
    return args

###############################################################################
# This function permit to get all txt file containned in input directory
###############################################################################


def get_genotype_files(input_directory):
    # Search for file ending with txt extension in a given directory
    input_files = glob.glob(input_directory + '*genotype*')
    return input_files

###############################################################################
# This function allow to determine for every results if it is a FP, TP or FN
###############################################################################


def result_analysis(result_file, pattern_size):
    result_type = 'FN'
    for line in result_file:
        if line[0]=='{':
            result_type = 'FP'
            pattern=line[line.find("{")+1:line.find("}")]
            if pattern.count('M') == pattern_size:
                result_type = 'TP'
                break
    return result_type

###############################################################################
# get the size of the pattern in simulated datas
###############################################################################


def get_pattern_size(input_files):
    file = open(input_files[0], 'r')
    size = file.readline().count('M')
    return size


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
    if (true_positive + false_negative)==0:
        recall = true_positive / 1
    else:
        recall = true_positive / (true_positive + false_negative)
    return recall

###############################################################################
# This function permit to calcul precision
###############################################################################


def calc_precision(true_positive, false_positive):
    if (true_positive + false_positive)==0:
        precision = true_positive / 1
    else:
        precision = true_positive / (true_positive + false_positive)
    return precision

###############################################################################
# This function permit to create measure file and calcul f_measure
###############################################################################


def calc_f_measure(recall, precision):
    f_measure = 2 / (1 + recall + 1 + precision)
    return f_measure

###############################################################################
# This function permit to powers file and calcul powers
###############################################################################


def calc_power(true_positive, number_of_execution):
    power = true_positive / number_of_execution
    return power

###############################################################################
# Main function
###############################################################################


def main():
    # Initialization of some variables
    scores = []
    # Get argument parser
    args = get_arguments()
    # Get argument passed to script
    input_directory = args.input
    output_directory = args.output
    number_of_execution = args.nruns
    method = args.method
    # Get list of genotype files in input directory
    input_files = get_genotype_files(input_directory)
    # Get pattern size from one file
    pattern_size = get_pattern_size(input_files)
    # Get number of files in input directory
    n_files = len(input_files)

    print(os.path.basename(input_directory))
    # For every files
    for file in input_files:
        pheno_file = file.replace("genotype", "phenotype")

        TP = 0
        FP = 0
        FN = 0

        for i in range(0, number_of_execution):
            if method == 'smmb_aco':
                os.system('./smmb_aco/smmb_aco.exe '+file+' '+pheno_file+' '+'./evaluation/parameters_smmb.txt')
            elif method == "vns":
                os.system('./vns/vns.exe '+file+' '+pheno_file+' '+'./evaluation/parameters_vns.txt')
            else:
                sys.exit("wrong method name, available methods : vns, smmb_aco")

            result_file = glob.glob("./evaluation/temp_results/" + "*")
            print(result_file)
            result_file1 = open(result_file[0], 'r')

            result_type = result_analysis(result_file1, pattern_size)
            result_file1.close()
            os.remove(result_file[0])
            if result_type == "TP":
                TP +=1
            elif result_type == "FP":
                FP +=1
            elif result_type == "FN":
                FN +=1

        recall = calc_recall(TP, FN)
        prec = calc_precision(TP, FP)
        f_measure = calc_f_measure(recall, prec)
        power = calc_power(TP, number_of_execution)
        scores.append((os.path.basename(os.path.normpath(file)), TP, FP, FN, recall, prec, f_measure, power))
        break
    power_file = open('./evaluation/result_eval/' + os.path.basename(os.path.normpath(input_directory)) +'_'+ method + ".csv", 'w')
    power_file.write("Filename,TP,FP,FN,recall,precision,f_measure,power"+"\n")
    for res in scores:
        power_file.write(str(res[0])+","+str(res[1])+","+str(res[2])+","+str(res[3])+","+str(res[4])+","+str(res[5])+","+str(res[6])+","+str(res[7])+"\n")






if __name__ == "__main__":
    main()
