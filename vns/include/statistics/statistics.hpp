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


    static void init_combinations();
    static void recursive_combination(unsigned step_val, unsigned array_index, unsigned _size_of_pattern, std::vector<unsigned> tuple, vector<unsigned> const& _possible_values, vector<vector<unsigned> > & _all_combinations);

    static void test();



private:

};

#endif
