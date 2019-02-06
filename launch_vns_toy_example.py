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
os.system('./vns/vns.exe ./vns/toy_example/toy_dataset/genotypes_toy_dataset.txt ./vns/toy_example/toy_dataset/phenotypes_toy_dataset.txt ./vns/toy_example/parameters.txt >> ./vns/log/log.log')
