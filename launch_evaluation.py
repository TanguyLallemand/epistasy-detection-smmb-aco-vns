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
os.system('./evaluation/eval_simu.py -i ./gametes_datas/datas/repository/model1/model1_0_01p_0005h_005m/ -o ./evaluation/result_eval/ -n 100 -m vns -nf 10 -l M')
os.system('./evaluation/eval_simu.py -i ./gametes_datas/datas/repository/model1/model1_0_01p_0005h_005m/ -o ./evaluation/result_eval/ -n 100 -m smmb -nf 10 -l M')
