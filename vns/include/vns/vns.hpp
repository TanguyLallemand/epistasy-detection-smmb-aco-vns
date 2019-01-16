/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#ifndef VNS_HPP
#define VNS_HPP

#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "global.hpp"
#include "statistics.hpp"

class vns
{

public:
    //==========================================================================
    // Constructor
    //==========================================================================
    vns(data_parsing dataset, parameters_parsing _params);
    //==========================================================================
    // Method
    //==========================================================================
    void run();
private:
    //==========================================================================
    // Variables
    //==========================================================================
    boost_matrix _genos_matrix;
    boost_vector_int _phenos_vector;
    boost_vector_string _snp_id;
    string _filename;
    string _output_directory;
    string _output_prefix;

    // Result map
    map<vector<unsigned>, vector<float>> _optimum_set;
    // Miscellaneous
    int _verbose;
    double _duration;
    // Parameters of algorithm
    unsigned _iteration_num;
    unsigned _pat_size_max;
    unsigned _pat_size_min;
    unsigned _k_max;
    unsigned _l_max;
    unsigned _max_it_vns;
    unsigned _max_it_local_search;
    float _alpha;
    // Seed for random
    std::mt19937 _rng;
    //==========================================================================
    // Methods
    //==========================================================================
    // Perturbation
    vector<unsigned> shake(vector<unsigned> pattern, unsigned k);
    // Local search
    vector<float> local_search(vector<unsigned> second_x, vector<unsigned> & third_x);
    void save_local_optimum(vector<unsigned> & x, vector<float> & x_score);
    // Output results in a file
    void write_result_file();
    vector<float> test_pattern(vector<unsigned> const& pattern);
    vector<unsigned> generate_starting_pattern();
    //==========================================================================
    // Miscellaneous
    //==========================================================================
    static bool compareFunc(pair<vector<unsigned>, vector<float>> const& a, pair<vector<unsigned>, vector<float>> const& b);
    void print_parameters();
};
#endif
