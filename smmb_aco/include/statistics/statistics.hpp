/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include "contingencies.hpp"
#include "global.hpp"

class statistics
{
public:
    //==========================================================================
    // Static functions to do chi 2 test and get associated score and p value
    //==========================================================================
    static float compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    static float compute_chi_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);
    //==========================================================================
    // Static functions to do chi 2 conditionnaly to variables test and get
    // associated score and p value
    //==========================================================================
    static boost_vector_float make_contingencies_chi_2_conditional_test_indep(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost_vector_int const& _phenos_vector, std::list<unsigned> const& cond_genos_indexes);
    static boost_vector_float compute_chi_2_conditional_test_indep(std::vector<contingencies> contingencies_vector, unsigned int liberty_degree);
    //==========================================================================
    // Tools for statitics
    //==========================================================================
    static unsigned compute_liberty_degree(boost_matrix_float const& contingency_table);
};

#endif
