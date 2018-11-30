# VNS documentation
## Variable Neightbor Search (VNS)
VNS is a program permitting to find epistasy patterns using Variable Neighbor Search.

## Compilation instructions
In order to compile this method a makefile was given at the root of VNS's folder. Please notice that a general makefile was generated at the project's root permitting to produce executables of both methods.
### Preparation of pre requist
In the makefile please change BOOST_FOLDER value with path of the installed BOOST library on current workstation.
Moreover, please check if current g++ version is compatible with C++11 functionalities (version >= 5.0). Please also check BOOST library. Authors cannot guaranty compatibility with other BOOST version even if BOOST version seems to be compatible. For this method version 1.61.0 was employed.
### Compilation
To compile this method, please call those lines a the root of the project:

    cd vns
    make

I a recompilation is needed use:

    cd smmb_aco
    make clean
    make clean

Execution of program with:

    ./vns.exe <path_to_genotype_dataset> <path_to_phenotype_dataset>

## Parameters
Parameters of this program are stored in parameters folder into a file called parameters.txt. This file is commented to permit an easy setup for particular needs. This allows to tweaks many parameters like output directory or output prefix and many parameters specific to VNS.
## Launch the analysis of the toy dataset
With this project is provided a python script allowing a simple launch on every platform equipped with python. This script will launch program on a naive dataset with default configuration given by authors.
To execute this script please call:

    ./launch_vns_toy_example.py

## Expected output
This method give as output a txt file gathering all results.
<!-- TODO mettre ce qu on attend -->
Those results can be checked using our evaluation tool given in evaluation folder.
Using toy dataset:
<!-- TODO mettre ce qu on attend -->

## Built With
-   [BOOST](https://www.boost.org/) - peer-reviewed portable C++ source libraries

## Authors
Tanguy Lallemand
Jonathan Cruard
