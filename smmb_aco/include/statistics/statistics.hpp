/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <vector>
#include <list>
#include <boost/math/distributions/chi_squared.hpp>

#include "global.hpp"

// liste des trucs à implémenter qui sont liés aux STATISTICS_HPP
    //test d'indépendance du chi 2
    //fonction pour trouver la pvalue d'une solution
    //test d'indépendance du chi 2 conditionnellement a une variable


class statistics
{
public:
    static float compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    static float compute_chi_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);
    static unsigned int compute_liberty_degree(boost_matrix_float const& contingency_table);
    static float statistics::make_contingencies_chi_2_conditional_test_indep(boost_matrix_float const& _genos_matrix, boost_vector_int const& _phenos_vector);
    static float statistics::compute_chi_2_conditional_test_indep();




private:

};

#endif
