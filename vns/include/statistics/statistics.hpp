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


    void get_all_combinations(boost_vector_int & subset, list<list<unsigned>> & combi_list);
    void generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset);
};

#endif
