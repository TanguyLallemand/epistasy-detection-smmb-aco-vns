#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
#
#
# Please call this script in a python environment with needed packages.
# In a folder called environment_to_execute_python a file is given to build a
# virtual environment containning all needed packages.
# Please read readme for more informations.
#
#
################################################################################
# Arguments                                                                    #
################################################################################
#
#
#-o OUTPUT, --output OUTPUT: Give an output directory
#-p PREFIX, --prefix PREFIX: Common prefix for files
#-f FILE, --file FILE  Number of file to generate
#-v VARIABLE, --variable VARIABLE: Number of variable to generate
#-pa CASE, --case CASE: Number of case to generate
#-c CONTROL, --control CONTROL: Number of control to generate
#-s SIZE_PATTERN, --size_pattern SIZE_PATTERN: Size of epistasis pattern to generate

import os
os.system('./simu_naive.py -p _simu_naive -f 1 -v 28 -pa 2000 -o toy_dataset_simu_naive -c 2000 -s 2')
