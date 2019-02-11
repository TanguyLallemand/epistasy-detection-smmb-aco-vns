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

################################################################################
# Creation of python environment                                               #
################################################################################

echo -e "\033[33m""Creating conda environment following configuration of workstation""\033[0m"
if which conda; then
    conda env create -f ./environment_to_execute_python/environment.yml
    else
    pip install virtualenv
    virtualenv -p /usr/bin/python3.6 projet_c
    source projet_c/bin/activate
    pip install -r environment.txt
    deactivate
fi
echo -e "\033[32m""Complete""\033[0m"
################################################################################
# Cleaning of compiled files, and recompile whole project                      #
################################################################################
echo -e "\033[33m""Cleaning compiled files""\033[0m"
make install
make clean
echo -e "\033[33m""Compilation of both methods""\033[0m"
make
echo -e "\033[33m""Project is ready to run, please use ./launch_vns_toy_example or ./launch_smmbaco_toy_example to use methods ""\033[0m"

echo -e "\033[33m""Do you want to generate a naive dataset using default parameters? ""\033[0m"
answer=("[y] yes" "[n]  no")
select choice in "${answer[@]}" ; do
    case $REPLY in
        1|y)
        echo -e "\033[33m""Runing a naive simulation""\033[0m"
        ./launch_simu_naive_toy_example.py
        echo -e "\033[32m""Complete""\033[0m"
        break
        ;;
        2|n)
        echo -e "\033[32m""Script will exit""\033[0m"
        break
        ;;
    esac
done
echo -e "\033[33m""Do you want to an evaluation on generated dataset using default parameters on both methods? ""\033[0m"
answer=("[y] yes" "[n]  no")
select choice in "${answer[@]}" ; do
    case $REPLY in
        1|y)
        echo -e "\033[33m""Runing an evaluation""\033[0m"
        ./evaluation/eval_simu.py -i ./toy_dataset_simu_naive/ -o ./evaluation/result_eval/ -n 50 -m vns -nf 50 -l X
        ./evaluation/eval_simu.py -i ./toy_dataset_simu_naive/ -o ./evaluation/result_eval/ -n 50 -m smmb_aco -nf 50 -l X
        echo -e "\033[32m""Complete""\033[0m"
        break
        ;;
        2|n)
        echo -e "\033[32m""Script will exit""\033[0m"
        break
        ;;
    esac
done
