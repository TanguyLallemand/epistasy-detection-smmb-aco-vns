#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
import os

# usage: eval_simu.py [-h] -i INPUT -o OUTPUT [-n NRUNS] -m METHOD
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Give an input directory
#   -o OUTPUT, --output OUTPUT
#                         Give an output directory
#   -n NRUNS, --nruns NRUNS
#                         Number of method executions to be performed on each
#                         file in the dataset
#   -m METHOD, --method METHOD
#                         Method to test

os.system('./evaluation/eval_simu.py -i ./smmb_aco/toy_example/toy_dataset/ -o ./evaluation/result_eval/ -n 100 -m smmb_aco')
