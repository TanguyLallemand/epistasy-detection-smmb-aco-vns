/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#ifndef PARAMETERS_PARSING_HPP
#define PARAMETERS_PARSING_HPP

#include "global.hpp"

class parameters_parsing
{
public:
    //==========================================================================
    // Constructor
    //==========================================================================
    parameters_parsing(string parameters_file);
    //==========================================================================
    // Method
    //==========================================================================
    void import_line(std::string const& line);
    //==========================================================================
    // Variables
    //==========================================================================
    // Name of files given in parameters
    string genos_file;
    string phenos_file;
    // Store asked output directory given in parameters
    string output_directory;
    // Store asked output file prefix given in parameters
    string output_prefix;
    int header;
    char separator;
    // Significance threshold
    float alpha;
    float precision;

    // Variables filled with parameters.txt
    //
    unsigned n_mbs;
    unsigned n_smmb_aco_runs;
    // Number of ACO iterations
    unsigned aco_n_iterations;
    // Number of ants
    unsigned aco_n_ants;
    // Size of the subset of sampled variables from _genos_matrix for each ants (so SNPs are sampled not individuals)
    unsigned aco_set_size;
    // Size of a combination of variables sampled from _subset_size
    unsigned subset_size_small;
    unsigned n_trials_to_learn_mbs;
    unsigned n_trials_to_learn_1_mb;
    // Initial pheromone value of each variable, at the beginning of the time perfect equality because we have no knowledge. To be treated as a vector
    float aco_tau_init;
    //Update of pheromone levels
    // Evaporation rate
    float aco_rho;
    // Values used in evaporation rates updates
    float aco_lambda;
    float aco_eta;
    // Two constants used to adjust the respective weights between pheromone levels and a priori knowledge. There may be SNPs that we know and therefore we give it a good rating. It therefore makes it possible to adjust the cursor between the importance of pheromones and the importance of knowledge
    float aco_alpha;
    float aco_beta;
private:
    vector<string> split(string const& s, char delim);
};

#endif
