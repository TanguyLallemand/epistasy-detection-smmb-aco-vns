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
    // Static functions to do g 2 test and get associated score and p value
    //==========================================================================
    static float compute_g_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);
    //==========================================================================
    // Static functions to do g 2 conditionnaly to variables test and get
    // associated score and p value
    //==========================================================================
    static boost_vector_float make_contingencies_g_2_conditional_test_indep(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost_vector_int const& _phenos_vector, std::list<unsigned> const& cond_genos_indexes);
    static boost_vector_float compute_g_2_conditional_test_indep(std::vector<contingencies> contingencies_vector, unsigned int liberty_degree);
};

#endif
