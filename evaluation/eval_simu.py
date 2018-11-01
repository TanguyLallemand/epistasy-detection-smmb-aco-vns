#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2 BB
#          Jonathan Cruard M2 BB
#
#
# Le script d'évaluation eval_simu.xxx (R ou Python) prendra en entrée :
# un nom de répertoire d'entrée
# un nom de répertoire de sortie
# le nombre d'exécutions de la méthode à réaliser sur chaque fichier du jeu de données simulées
# (nommé n_runs ci après)
# Pour chaque fichier <identifiant_fichier_i.txt> du jeu de données simulées, un fichier
# <identifiant_fichier_i>_results.txt de n_runs mots dans {TP, FP, FN} sera généré.


# generation de données simulées via script.
# Sortie:
#     SNP et score
#     {var1, var2, var3} <score>
#     {var1, var45, var2000, var5000} <score>
#     ...
#     {var67, var340} <score>
#
# On doit pouvoir determine Faux positifs et faux negatifs afin de pouvoir determiner recall precision.
#
# Exemple de fichier <identifiant_fichier_i>_results.txt :
# TP
# FN
# TP
# FP
# ...
# FN
#
#
# Pour chaque jeu de données (comportant par exemple n_files fichiers), un fichier
# f_measures.txt sera généré. Il comportera les n_files f-measures calculées à partir des n_runs fichiers <identifiant_fichier_i>_results.txt générés pour chacun # des fichiers <identifiant_fichier_i.txt> (1 ≤ i ≤ n_files) :
# recall = #TP ⁄ (#TP + #FN)
# precision = #TP ⁄ (#TP + #FP).
# f-measure =2 / (1 ⁄ recall +1 ⁄ precision)
# Pour chaque jeu de données (comportant par exemple n_files fichiers), un fichier powers.txt
# sera généré. Il comportera les n_files f-measures calculées à partir des n_runs fichiers
# <identifiant_fichier_i>_results.txt générés pour chacun des fichiers <identifiant_fichier_i.txt>
# (1 ≤ i ≤ n_files).
# power = #TP / n_run
#
# Regle pour pattern de taille 2
# si contient pattern simulé ca va donc TP true Positive
#
# Faux positif (on a trouve un truc faux)
# Faux negatif (on a pas trouve)
#     Fichier vide
#     Pas la bonne taille (donc la inferieur a 2)
#
# pour avoir un seul resultat par fichier, si on trouve un TP tout le fichier est TP
# sinon c'est a la majorité
#
# f measure pour eviter qu il plante a cause division par 0, si y a un seul TP et le reste en FP ou FP ca passe mais si pas de TP et tt en FN alors division par 0. Si TP=0 alors pb. Dans ce cas la j en enleverai un a celui qui peut pas etre calcule et je le met a l autre qui est a 0

import argparse
import glob
import os
import numpy as np


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


def get_input_files(input_directory):
    # Search for file ending with txt extension in a given directory
    input_files = glob.glob(input_directory + '*.txt')
    return input_files


def parsing_result_file(result_file):
    # Parsing of result file
    result = []
    for line in result_file:
        if line == pattern:
            result.append('TP')
            true_positive += 1
        elif line.strip() or len(line) != len(pattern):
            result.append('FN')
            false_negative += 1
        elif line != pattern:
            result.append('FP')
            false_positive += 1
    return [result, true_positive, false_negative, false_positive]


def creation_of_output_file(output_directory, results_file_parsed):
    base = os.path.basename(name_of_file)
    name = os.path.splitext(base)[0] + '_results.txt'
    if os.path.isdir(output_directory):
        with open(name, 'w') as output_file:
            pass


def creation_of_measure_file(recall, precision):
    with open('f_measures.txt', 'w') as measure_file:
        f_measure = 2 / (1 ⁄ recall + 1 ⁄ precision)
        measure_file.write(f_measure)



def creation_of_powers_file():
    with open('powers.txt', 'w') as powers_file:


        # Get argument parser
args = get_arguments()
input_directory = args.input
output_directory = args.output
number_of_execution = args.nruns

input_files = get_input_files(input_directory)
n_files = len(input_files)

for file in input_files:
    with open(file, 'r') as result_file:
        results_file_parsed = parsing_result_file(result_file)
        array_of_results_string = results_file_parsed[0]
        true_positive = results_file_parsed[1]
        false_negative = results_file_parsed[2]
        false_positive = results_file_parsed[3]
        creation_of_output_file(output_directory, array_of_results_string)
        recall = calc_recall(true_positive, )
        precision = calc_precision(true_positive, )

        creation_of_measure_file(recall, precision)
    # creation_of_powers_file():
