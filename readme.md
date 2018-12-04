# Detection of epistasis patterns in genetic data

## Methods used for this project

In this project two methods were implemented to search for epistasy patterns in genetic data.

-   SMMB-ACO for Multiple Markov Blankets by Ant Colony Optimization.
-   VNS for Variable Neighborhood Search.

## Compilation instructions

In order to compile the whole project a makefile is provided at project's root. This makefile will call each method's makefile and produce both executables.

### Preparation of prerequisite

In makefiles please change BOOST_FOLDER value with path of the installed boost library on current workstation.
Moreover, please check if current g++ version is compatible with C++11 functionalities (version >= 5.0). Please also check boost library's version. Authors cannot guaranty compatibility with other boost version even if boost version seems to be compatible. For this method version 1.61.0 was employed.

### Compilation

To compile this method, please call this lines at the root of the project:

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

It is possible to call this tool using a default call using launch_simu_naive_toy_example.py. In this launcher is also added explications about arguments used in call. It is also possible to display help using

    simu_naive.py --help

## Expected output

On each data set constituted by one genotype file and his associated phenotype file epistasis pattern is searched using both methods. Here is an example of expected output file.

    # Result from vns
    # Pattern || occurences || chi2-score || p-value || unreliable case
    {N2,N12,M0P7} || 1 || 43.0154 || 0.0192421 || 12
    {N7,N20,M0P8} || 1 || 48.6441 || 0.00454652 || 9
    {N8,N11,N24} || 1 || 44.8019 || 0.0123739 || 10
    {N11,N24} || 1 || 17.463 || 0.0256333 || 0
    {N12,N24,M0P7} || 1 || 42.5037 || 0.0217698 || 4
    {N13,M0P7,M0P8} || 1 || 99.0419 || 0 || 3
    # Time of execution: 37 seconds

This method give as output a text file gathering all results. Those results can be checked using our evaluation tool given in evaluation folder.
Using a toy data set:

<!-- TODO add  -->

## Packages and library used in this project

-   [BOOST](https://www.boost.org/) a peer reviewed c++ collection of library.
-   [NumpY](http://www.numpy.org/) a python library to handle with matrix.
-   [Pandas](https://pandas.pydata.org/) a Python Data Analysis Library.

## Authors

Tanguy Lallemand -
Jonathan Cruard
