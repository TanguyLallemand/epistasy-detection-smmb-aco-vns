/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include "contingencies.hpp"
#include "global.hpp"

class statistics
{
public:


    //==========================================================================
    // Static functions to do g2 test and get associated score and p value
    //==========================================================================


    static float compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    static float compute_g2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);


    //==========================================================================
    // Tools for statitics
    //==========================================================================


    static unsigned compute_liberty_degree(boost_matrix_float const& contingency_table);


    //==========================================================================
    // Tools for generate contingencies tables
    //==========================================================================


    static void init_combinations(unsigned size_of_pattern, vector<vector<unsigned> > & all_combinations, vector<unsigned> & possible_values);
    static void recursive_combination(unsigned step_val, unsigned array_index, std::vector<unsigned> tuple, unsigned size_of_pattern, vector<vector<unsigned> > & all_combinations, vector<unsigned> & possible_values);

};

#endif
