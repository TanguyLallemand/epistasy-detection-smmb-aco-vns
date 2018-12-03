# Detection of epistasis patterns in genetic data
## Methods used for this project
In this project two methods were implemented to search for epistasy patterns in genetic data.

- SMMB-ACO for Multiple Markov Blankets by Ant Colony Optimization.
- VNS for Variable neighborhood search.
## Compilation instructions
In order to compile the whole project a makefile is provided at project's root. This makefile will call each of the makefile methods to build a program containing the two implemented methods.
### Preparation of pre requist
In makefile please change BOOST_FOLDER value with path of the installed boost library on current workstation.
Moreover, please check if current g++ version is compatible with C++11 functionalities (version >= 5.0). Please also check boost library. Authors cannot guaranty compatibility with other boost version even if boost version seems to be compatible. For this method version 1.61.0 was employed.
### Compilation
To compile this method, please call this lines at the root of the project:


    make


If a recompilation is needed please use:

    make clean

Execute program with:

<!-- TODO mettre ici le call de smmb -->
## Parameters
It is possible to tweaks differents parameters for both methods. To change parameters of SMMB ACO method use parameters.txt localised in smmb_aco/paramters. This file is commented to permit an easy setup for particular needs. If some tweaks are needed for VNS method,
<!-- TODO mettre le chemin des parametres -->
## Launch the analysis of the toy dataset
With this project is provided a python script allowing a simple launch on every platform equiped with python. This script will launch program on a naive dataset with default configuration given by authors.
To execute this script please call:

    ./launch_toy_dataset_analysis.py

## Expected output
This method give as output a txt file gathering all results. Those results can be checked using our evaluaton tool given in evaluation folder.
Using toy dataset:
<!-- TODO mettre ce qu on attend -->
## Packages and library used in this project
- [BOOST](https://www.boost.org/) a peer reviewed c++ collection of library.
- [NumpY](http://www.numpy.org/) a python library to handle with matrix.

## Authors
Tanguy Lallemand -
Jonathan Cruard
