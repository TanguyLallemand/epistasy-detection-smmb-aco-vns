#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
#
#
################################################################################
# Arguments                                                                    #
################################################################################
#
#
# usage: eval_simu.py [-h] -i INPUT -o OUTPUT [-l CAUSAL_ID] -m METHOD
#                     [-n NRUNS] [-nf NFILES]
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Give an input directory
#   -o OUTPUT, --output OUTPUT
#                         Give an output directory
#   -l CAUSAL_ID, --causal_id CAUSAL_ID
#                         Letter to identify causal SNPs
#   -m METHOD, --method METHOD
#                         Method to test
#   -n NRUNS, --nruns NRUNS
#                         Number of method executions to be performed on each
#                         file in the dataset
#   -nf NFILES, --nfiles NFILES
#                         Number of file to test in the dataset

import os
os.system('./evaluation/eval_simu.py -i ./simu_naive_CARLUER_OUEDRAOGO/simu_naive_2snp_0_059_0_25/ -o ./evaluation/result_eval/ -n 100 -m vns -nf 100 -l M')
os.system('./evaluation/eval_simu.py -i ./simu_naive_CARLUER_OUEDRAOGO/simu_naive_2snp_0_059_0_25/ -o ./evaluation/result_eval/ -n 100 -m smmb_aco -nf 100 -l M')

os.system('./evaluation/eval_simu.py -i ./Simu_naive/Simu_naive_2snp_0.25/ -o ./evaluation/result_eval/ -n 100 -m vns -nf 100 -l X')
os.system('./evaluation/eval_simu.py -i ./Simu_naive/Simu_naive_2snp_0.25/ -o ./evaluation/result_eval/ -n 100 -m smmb_aco -nf 100 -l X')
