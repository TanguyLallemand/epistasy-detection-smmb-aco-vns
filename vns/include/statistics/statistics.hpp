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
    // Static functions to do g2 test and get associated score, p value and
    // Number of non reliable tests
    //==========================================================================
    static  vector<float> compute_p_value(vector<boost::numeric::ublas::matrix_column<boost_matrix>> const& pattern_datas, boost_vector_int const& _phenos_vector);
    static  vector<float> compute_g2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);
    //==========================================================================
    // Tools for statitics, compute liberty degree following contingency table's
    // Size
    //==========================================================================
    static unsigned compute_liberty_degree(boost_matrix_float const& contingency_table);
    //==========================================================================
    // Tools to generate contingencies tables in recursive fashion
    //==========================================================================
    static vector<vector<unsigned>> init_combinations(unsigned pattern_size);
    static void recursive_combination(unsigned _size_of_pattern, std::vector<unsigned> tuple, vector<unsigned> const& _possible_values, vector<vector<unsigned>> & _all_combinations);
};
#endif
