#! /usr/bin/python3
#-*-coding: utf-8-*-
# Louison Fresnais M2BB
# François Courtin M2BB

#####################################################################
#   Cette fonction permet le parsing des paramètres du script		#
#####################################################################

def argument_parsing():

    import argparse 

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Naive data simulation script",
        epilog="""The objective of this script is to create naively simulated data for an
            epistasia pattern detection algorithm.
            The output of this script is:  
                - a number of simulated genotype data
                - a number of simulated phenotype data
            
            Input parameters are the following: 
            -1: Number of files to generate per run
            -2: A run prefix to identify each result file from this run 
            -3: The number of variables to simulate
            -4: The number of case to simulate
            -5: The number of controls to simulate
            Missing parameters will be asked in the console.

            Contact: fresnaislouison@gmail.com""")

    parser.add_argument('-nf','--file', type=int, help='enter the number of file you want to generate in one run')
    parser.add_argument('-pf','--prefix', type=str, help='enter a disctinctive prefix that will be applied to each result file for this run')
    parser.add_argument('-nv','--variable', type=int, help='enter the number of variables to simulate')
    parser.add_argument('-nc','--case', type=int, help='enter the number of case to simulate')
    parser.add_argument('-no','--control', type=int, help='enter the number of controls to simulate')
    parser.add_argument('-pa','--pattern', type=int, help='enter the size of the epistasia pattern')

    args = parser.parse_args()

    return args

