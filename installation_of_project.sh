#!/bin/bash

################################################################################
# Creation of directories                                                      #
################################################################################
echo -e "\033[33m""Creation of architecture""\033[0m"
# Check for architecture
mkdir -p ./evaluation/
mkdir -p ./evaluation/result_eval
mkdir -p ./evaluation/temp_results
mkdir -p ./smmb_aco/obj
mkdir -p ./vns/obj
echo -e "\033[32m""Complete""\033[0m"
echo -e "\033[33m""Creating conda environment""\033[0m"
conda env create -f ./environment_to_execute_python/environment.yml
echo -e "\033[32m""Complete""\033[0m"
echo -e "\033[33m""Cleaning compiled files""\033[0m"
make install
make clean
echo -e "\033[33m""Compilation of both methods""\033[0m"
make
echo -e "\033[33m""Project is ready to run, please use ./launch_vns_toy_example or ./launch_smmbaco_toy_example ""\033[0m"
