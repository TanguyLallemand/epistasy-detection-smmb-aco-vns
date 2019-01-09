/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP

#include "global.hpp"

#include "file_parsing.hpp"
#include "statistics.hpp"
#include "tools.hpp"
#include "parameters_parsing.hpp"

class smmb_aco
{
public:
    //==========================================================================
    // Constructor
    //==========================================================================
    smmb_aco(data_parsing dataset, parameters_parsing _params);
    // running SMMB_ACO
    void run();

private:
    //==========================================================================
    // Get object given as parameters
    //==========================================================================
    boost_matrix _genos_matrix;
    boost_vector_int _pheno_vector;
    boost_vector_string _snp_id;
    string _filename;
    // Miscellaneous
    int _verbose;
    //==========================================================================
    // Variables initialized by the constructor from parameters
    //==========================================================================
    //variable to now if a second pass is asked and on which pass we are
    unsigned _pass_number;
    // Number of ACO iterations
    unsigned _n_it_n;
    // Number of ants
    unsigned _n_ant;
    // Size of the subset of sampled variables from _genos_matrix for each ants (so SNPs are sampled not individuals)
    unsigned _subset_size;
    // Size of a combination of variables sampled from _subset_size
    unsigned _sub_subset_size;
    // Maximum number of iterations to explore the research space
    unsigned _n_it; //TODO check tout ça car tu a 3 trucs de limite dans les paramètres et on en utilise que 2
    // Significance threshold
    float _alpha_stat;
    // Initial pheromone value of each variable, at the beginning of the time perfect equality because we have no knowledge. To be treated as a vector
    float _tau_0;
    //Update of pheromone levels
    // Evaporation rate
    float _rho;
    // Values used in evaporation rates updates
    float _lambda;
    // Two constants used to adjust the respective weights between pheromone levels and a priori knowledge. There may be SNPs that we know and therefore we give it a good rating. It therefore makes it possible to adjust the cursor between the importance of pheromones and the importance of knowledge
    float _alpha_phero;
    float _beta_phero;
    // Store asked output directory given in parameters
    string _output_directory;
    // Store asked output file prefix given in parameters
    string _output_prefix;
    // Init a variable to hold time of execution
    double _duration;
    //rng seed
    std::mt19937 _rng;

    //==========================================================================
    // Structure used in this algorithm
    //==========================================================================
    // Global memory of processed tests through 1 iteration
    std::map<unsigned, list<float>> _mem;
    // Every case of the vector is the memory of 1 ant
    boost::numeric::ublas::vector<map<unsigned, list<float>>> _mem_ant;
    // 1 case is best MB found by the corresponding ant
    boost::numeric::ublas::vector<list<unsigned>> _markov_blanket_a;
    // Stocking the result of each run with pattern as key and occurence count of it
    map<list<unsigned>, unsigned> _markov_blanket_s;
    // Associated vector of scores
    boost::numeric::ublas::vector<boost_vector_float> _stats_results;
    //Vectors about pheromons
    boost_vector_float _eta;
    boost_vector_float _tau;
    boost_vector_float _pheromone_distrib;
    //==========================================================================
    // Main functions used in this algorithm
    //==========================================================================

    void learn_MB(boost_vector_int & ant_subset, list<unsigned> & MB_a_ref, std::map<unsigned, list<float>> & mem_ant_ref);
    void forward(list<unsigned> & MB_a_ref, boost_vector_int const& ant_subset, std::map<unsigned, list<float>> & mem_ant_ref);
    void backward(list<unsigned> & MB_a_ref);

    //==========================================================================
    //Functions for modularity
    //==========================================================================
    // Permit to add pheromones on a good SNP
    void update_tau();
    // Compute sub_subset
    void sub_sampling(boost_vector_int & sub_subset, boost_vector_int const& ant_subset, list<unsigned> MB_a_ref);
    // Get all combinations of a given subset
    void get_all_combinations(boost_vector_int & sub_subset, list<list<unsigned>> & combi_list);
    // Generate all combinations of a subset in a recursive fashion
    void generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset);
    // Determine best combination of SNP using g 2 conditional test of independance
    boost_vector_float best_combination(list<unsigned> & best_pattern, list<list<unsigned>> const& pattern_list, list<unsigned> & MB_a_ref, std::map<unsigned, list<float>> & mem_ant_ref);
    // Update pheromons using _tau and _eta
    void update_pheromon_distrib();
    // prepare next pass and lunch it
    void next_pass(list<unsigned> new_set);
    // Compute g 2 conditional test of independance on markov blanket
    void score_for_final_results();
    // Permit to show current results, usefull in debug process
    void show_results();
    // save the current iteration results into global variables
    void save_iteration_result();
    // Write results in a file
    void save_results();



    static bool compareFunc(pair<unsigned, float> const& a, pair<unsigned, float> const& b);
};
#endif
