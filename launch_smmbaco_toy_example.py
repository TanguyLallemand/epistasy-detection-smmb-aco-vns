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
# First: path of genotype data set
# Second: path of associated phenotype data set
# Third: path to parameters file
# Finally > is used to redirect standard output in a log file
import os
os.system('./smmb_aco/smmb_aco.exe ./smmb_aco/toy_example/toy_dataset/genotypes_toy_dataset.txt ./smmb_aco/toy_example/toy_dataset/phenotypes_toy_dataset.txt ./smmb_aco/toy_example/parameters.txt >> ./smmb_aco/outputs/log.log')
