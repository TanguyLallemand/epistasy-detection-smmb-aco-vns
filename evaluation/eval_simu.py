#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2 BB
#          Jonathan Cruard M2 BB


import glob
import numpy as np
import os
import time

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
        "-l", "--causal_id", help="Letter to identify causal SNPs", type=str, action='store')
    parser.add_argument(
        "-m", "--method", help="Method to test", type=str, action='store', required=True)
    parser.add_argument(
        "-n", "--nruns", help="Number of method executions to be performed on each file in the dataset", type=int)
    parser.add_argument(
        "-nf", "--nfiles", help="Number of file to test in the dataset", type=int)
    args = parser.parse_args()
    return args

###############################################################################
# This function permit to get all txt file containned in input directory
###############################################################################


def get_genotype_files(input_directory):
    # Search for file ending with txt extension in a given directory
    input_files = glob.glob(input_directory + '*geno*')
    if len(input_files)<1:
        input_files = glob.glob(input_directory + '*Geno*')
    return input_files

###############################################################################
# This function allow to determine for every results if it is a FP, TP or FN
###############################################################################


def result_analysis(result_file, pattern_size, causal_id):
    result_type = 'FN'
    for line in result_file:
        if line[0]=='{':
            result_type = 'FP'
            pattern=line[line.find("{")+1:line.find("}")]
            if pattern.count(causal_id) == pattern_size:
                result_type = 'TP'
                break
    return result_type

###############################################################################
# get the size of the pattern in simulated datas
###############################################################################


def get_pattern_size(input_files, causal_id):
    file = open(input_files[0], 'r')
    size = file.readline().count(causal_id)
    return size


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
# functions to write the different result files
###############################################################################
def write_f_measure(scores, input_directory, method):
    power_file = open('./evaluation/result_eval/' + os.path.basename(os.path.normpath(input_directory)) +'_'+ method + "_f_measure.csv", 'w')
    power_file.write("Filename,f_measure"+"\n")
    for res in scores:
        power_file.write(str(res[0])+","+str(res[6])+"\n")

def write_power(scores, input_directory, method):
    power_file = open('./evaluation/result_eval/' + os.path.basename(os.path.normpath(input_directory)) +'_'+ method + "_power.csv", 'w')
    power_file.write("Filename,power"+"\n")
    for res in scores:
        power_file.write(str(res[0])+","+str(res[7])+"\n")

def write_time(scores, input_directory, method):
    power_file = open('./evaluation/result_eval/' + os.path.basename(os.path.normpath(input_directory)) +'_'+ method + "_time.csv", 'w')
    power_file.write("Filename,average_time_per_run"+"\n")
    for res in scores:
        power_file.write(str(res[0])+","+str(res[8])+"\n")

def write_global_results(scores, start, end, input_directory, method):
    power_file = open('./evaluation/result_eval/' + os.path.basename(os.path.normpath(input_directory)) +'_'+ method + ".csv", 'w')
    power_file.write("Filename,TP,FP,FN,recall,precision,f_measure,power,average_time_per_run"+"\n")
    for res in scores:
        power_file.write(str(res[0])+","+str(res[1])+","+str(res[2])+","+str(res[3])+","+str(res[4])+","+str(res[5])+","+str(res[6])+","+str(res[7])+","+str(res[8])+"\n")
    power_file.write('### Total evaluation time : ' + str(end-start))

###############################################################################
# Main function
###############################################################################


def main():
    start = time.time()
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
    pattern_size = get_pattern_size(input_files, args.causal_id)
    # Get number of files in input directory
    n_files = len(input_files)

    # Global counter
    z = 0
    zmax = number_of_execution * min(args.nfiles, n_files)
    # For every files
    for file in input_files:
        start_file = time.time()
        pheno_file = file.replace("geno", "pheno")
        pheno_file = pheno_file.replace("Geno", "Pheno")
        TP = 0
        FP = 0
        FN = 0

        for i in range(1, number_of_execution):
            z+=1
            print(file + " Iteration " + str(i) + ' ' + str(z) + '/' + str(zmax))
            if method == 'smmb_aco':
                os.system('./smmb_aco/smmb_aco.exe '+file+' '+pheno_file+' '+'./evaluation/parameters_smmb.txt >> ./evaluation/smmb_aco_log.log')
            elif method == "vns":
                os.system('./vns/vns.exe '+file+' '+pheno_file+' '+'./evaluation/parameters_vns.txt >> ./evaluation/vns_log.log')
            else:
                sys.exit("wrong method name, available methods : vns, smmb_aco")

            result_file_name = glob.glob("./evaluation/temp_results/" + "*")
            result_file_handler = open(result_file_name[0], 'r')

            result_type = result_analysis(result_file_handler, pattern_size, args.causal_id)
            result_file_handler.close()
            os.remove(result_file_name[0])
            if result_type == "TP":
                TP +=1
            elif result_type == "FP":
                FP +=1
            elif result_type == "FN":
                FN +=1
        end_file = time.time()
        recall = calc_recall(TP, FN)
        prec = calc_precision(TP, FP)
        f_measure = calc_f_measure(recall, prec)
        power = calc_power(TP, number_of_execution)
        scores.append((os.path.basename(os.path.normpath(file)), TP, FP, FN, recall, prec, f_measure, power, str((end_file-start_file)/number_of_execution)))
        if args.nfiles == 0:
            break
        else:
            args.nfiles -= 1

    end = time.time()
    write_f_measure(scores, input_directory, method)
    write_power(scores, input_directory, method)
    write_time(scores, input_directory, method)
    write_global_results(scores, start, end, input_directory, method)



if __name__ == "__main__":
    main()
