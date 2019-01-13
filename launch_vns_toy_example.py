#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Tanguy Lallemand M2BB
#          Jonathan Cruard M2BB
import os

# This programm require three parameters:
# First: path of genotype data set
# Second: path of associated phenotype data set
# Third: path to parameters file
os.system('./vns/vns.exe ./gametes_datas/datas/repository/model1/model1_0_01p_0005h_005m/genotype_model1_0_01p_0005h_005m_005.csv ./gametes_datas/datas/repository/model1/model1_0_01p_0005h_005m/phenotype_model1_0_01p_0005h_005m_005.csv ./vns/toy_example/parameters.txt')
