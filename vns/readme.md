# VNS documentation
## Variable Neightbor Search (VNS)
VNS is a program permitting to find epistasy patterns using Variable Neighbor Search.

## Compilation instructions
In order to compile this method a makefile was given at the root of VNS's folder. Please notice that a general makefile was generated at the project's root permitting to produce executables of both methods.
### Preparation of pre requist
In the makefile please change BOOST_FOLDER value with path of the installed BOOST library on current workstation.
Moreover, please check if current g++ version is compatible with C++11 functionalities (version >= 5.0). Please also check BOOST library. Authors cannot guaranty compatibility with other BOOST version even if BOOST version seems to be compatible. For this method version 1.61.0 was employed.
### Compilation
To compile this method, please call those lines at the root of the project:

    cd vns
    make

I a recompilation is needed use:

    cd smmb_aco
    make clean
    make

Execution of program with:

    ./vns.exe <path_to_genotype_dataset> <path_to_phenotype_dataset> <path_to_parameters.txt>

## Summary of files of this project
| Module         | Files              | Usability                                                                                                  | Inputs                                                                                                                                                         | Outputs                                                                                       |
|----------------|--------------------|------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| Parsing        | file_parsing       | Allows to parse file containing data set                                                                   | Files names given as argument in program's call                                                                                                                |  One boost matrix containing all genotype data <br> One boost vector containing all phenotype data <br> boost vector with SNPs IDs |
| Parsing        | parameters_parsing | Allows to parse parameters and save them in parameters class                                               | Will use ./parameters/parameters.txt                                                                                                                           | Class object with all parameters as class variables                                           |
| Smmb aco       | smmb_aco           | Core of this project allows to run smmb-aco's algorithm                                                    | Parameters object containing all parameters <br> Data set object containing: <br>  - Matrix of genotype data <br>   - Vector of phenotype data               | Final markov blanket                                                                          |
| Statistics     | statistics         | Used to do g2 conditional test of independence                                                            | Matrix column of genotype (subset of genotype matrix with only tested SNPs)<br> Vector of phenotype Indexes of tested SNPs                                         |  g2 score and associated p-value Number of cell considered as non reliable because n<5         |
| Statistics     | contingencies      |  Build observed and expected contingencies table,   some functions useful to work on contingencies tables | Matrix column of genotype (subset of genotype matrix with only tested SNPs) <br> Vector of phenotype <br>  Indexes of tested SNPs Number of observation in tested subset |  One observed contingency table <br> One expected contingency table                                |
| Output writing | output_writing     | Can save algorithm's results in file                                                                       | Markov blanket and associated score                                                                                                                            | Final result file                                                                             |
| Tools          | tools              |   Used to generate random discrete distribution used to pick   SNPs when building subset                   |  Weight vector size of subset in construction <br> Random seed                                                                                                      | Subset of SNPs for every ant                                                                  |
## Parameters
Parameters of this program are stored in parameters folder into a file called parameters.txt. This file is commented to permit an easy setup for particular needs. This allows to tweaks many parameters like output directory or output prefix and many parameters specific to VNS.
## Launch the analysis of the toy data set
With this project is provided a python script allowing a simple launch on every platform equipped with python. This script will launch program on a naive data set with default configuration given by authors.
To execute this script please call:

    ./launch_vns_toy_example.py

## Expected output
This method give as output a txt file gathering all results.
    # Result from vns 
    # Pattern || occurences || chi2-score || p-value || unreliable case
    {N0,N18} || 1 || 19.8798 || 0.0108003 || 0
    {N3,N21,M0P7} || 1 || 52.6552 || 0.00149465 || 14
    {N4,N17,M0P8} || 1 || 49.0533 || 0.00407094 || 4
    {N6,N9} || 1 || 18.9766 || 0.0149854 || 0
    {N6,N21,M0P7} || 1 || 39.7037 || 0.0417014 || 14
    {N9,N10,N17} || 1 || 41.8464 || 0.0254582 || 10
    {N9,N20} || 1 || 19.9661 || 0.0104652 || 3
    {N10,M0P7} || 1 || 26.5278 || 0.000852644 || 0
    {N11,M0P7} || 1 || 29.0766 || 0.000307441 || 0
    {N12,M0P8} || 1 || 29.5862 || 0.000250101 || 0
    {N20,N21,M0P7} || 1 || 42.3177 || 0.022761 || 10
    {N21,M0P7} || 1 || 26.6079 || 0.000826061 || 4
    # Execution time : 13 seconds
Those results can be checked using our evaluation tool given in evaluation folder.
Using a toy data set:

Filename,TP,FP,FN,recall,precision,f_measure,power
genotypes_toy_dataset.txt,2,55,43,0.044444444444444446,0.03508771929824561,0.9617547806524185,0.02

## Built With
-   [BOOST](https://www.boost.org/) - peer-reviewed portable C++ source libraries

## Authors
Tanguy Lallemand -
Jonathan Cruard
