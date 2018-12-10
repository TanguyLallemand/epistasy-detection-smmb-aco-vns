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
    // Static functions to make a vector of contingency table and compute
    // associated theorical table
    //==========================================================================
    static boost_vector_float make_contingencies_g_2_conditional_test_indep(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost_vector_int const& _phenos_vector, std::list<unsigned> const& cond_genos_indexes);
    //==========================================================================
    // Static functions to compute g2 test conditionnaly to variables and get
    // associated score, p value and number of non reliable tests
    //==========================================================================
    static boost_vector_float compute_g_2_conditional_test_indep(std::vector<contingencies> contingencies_vector, unsigned int liberty_degree);
    //==========================================================================
    // Static functions to compute g2 test and get associated score
    //==========================================================================
    static float compute_g_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);
};

#endif
