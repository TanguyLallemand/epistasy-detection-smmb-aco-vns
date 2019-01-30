# Detection of epistasis patterns in genetic data

## Methods used for this project

In this project two methods were implemented to search for epistasy patterns in genetic data.

-   SMMB-ACO for Multiple Markov Blankets by Ant Colony Optimization.
-   VNS for Variable Neighborhood Search.

## Compilation instructions

In order to compile the whole project a makefile is provided at project's root. This makefile will call each method's makefile and produce both executables.

### Preparation of prerequisite

In makefiles please change BOOST_FOLDER value with path of the installed boost library on current workstation.
Moreover, please check if current g++ version is compatible with C++11 functionalities (version >= 5.0). Please also check boost library's version. Authors cannot guaranty compatibility with other boost version even if boost version seems to be compatible between them. For this method version 1.61.0 was employed.

### Compilation

To compile both methods, please call this lines at the root of the project:

    make

If a recompilation is needed please use:

    make clean

## Parameters

It is possible to tweaks different parameters for both methods. To change parameters of SMMB ACO method use [parameters.txt](smmb_aco/parameters/parameters.txt) localized in smmb_aco/parameters. This file is commented allowing an easy setup for particular needs. If some tweaks are needed for VNS method, please modify vns/parameters/parameters.txt. Please note that a particular section has been included in associated project report about some parameters and their effects.

## Generate a naive data set

This project is provided with a tool to generate some naive data set. This tool follows a logistic regression and produce a data set of given sizes with a given number of causal SNPs. This operation is performed by simu_naive.py.

### Prerequisite, generation of right virtual environment

Some packages are required to run this tool. List of them is provided with version in environment.yml file. In order to generate a virtual environment where simu_naive.py can be executed please follow next steps.

On a workstation with conda installed please execute this line to generate an environment called projet_c:

    conda env create -f environment.yml

Following this, activate this environment and execute simu_naive with arguments:

    source activate projet_c
    ./simu_naive.py -p simu_naive -f 1 -v 28 -pa 2000 -o toy_dataset -c 2000 -s 2
    source deactivate

### Run simu_naive.py

It is possible to call this tool using a default call using

    ./launch_simu_naive_toy_example.py

In this launcher is also added explications about arguments used in call. It is also possible to display help using

    simu_naive.py --help


## Expected output

On each data set constituted by one genotype file and his associated phenotype file epistasis pattern is searched using both methods. Here is an example of expected output file :

    # Result from vns
    # Pattern || occurences || chi2-score || p-value || unreliable case
    {SNP-N-1,SNP-N-17,SNP-C-0} || 1 || 422.528 || 0 || 0
    {SNP-N-3} || 1 || 6.7944 || 0.0334668 || 0
    {SNP-N-6,SNP-C-0} || 1 || 393.468 || 0 || 0
    {SNP-N-8,SNP-C-1} || 1 || 3898.56 || 0 || 0
    {SNP-N-19,SNP-N-20,SNP-C-0} || 1 || 413.348 || 0 || 0
    {SNP-C-0} || 1 || 389.899 || 0 || 0
    {SNP-C-1} || 2 || 3893.06 || 0 || 0
    # Time of execution:0seconds

Each method can be evaluated using the eval_simu.py script, it uses the method as many times as it is asked on each data file in the directory given as an argument. A file summarizing the results is produced, it includes the counting of the type of results (true positive: TP, false positive: FP, false negative: FN) as well as the precision calculations, f-measurement, etc.

#### Run eval_simu.py

This script allows to execute a given number of times the methods on files contained in a given folder. In order to execute this script with a default configuration please call it with following line:

    ./launch_evaluation.py

This launcher will perform each methods a given number of times on some files generated using gametes. If a particular configuration is needed, please check call used in launcher or display help using:

    ./evaluation/eval_simu.py --help

It will produce a file of this type :

    Filename,TP,FP,FN,recall,precision,f_measure,power
    genotypes_toy_dataset.txt,2,55,43,0.044444444444444446,0.03508771929824561,0.9617547806524185,0.02


## Packages and library used in this project

-   [BOOST](https://www.boost.org/) a peer reviewed c++ collection of library.
-   [NumpY](http://www.numpy.org/) a python library to handle with matrix.
-   [Pandas](https://pandas.pydata.org/) a Python Data Analysis Library.

## Authors

Tanguy Lallemand -
Jonathan Cruard
