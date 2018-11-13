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
    // constructor
    statistics(boost_matrix _genos_matrix, boost_vector_int _phenos_vector, parameters_parsing _params);
    static boost_matrix_float make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    static float compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    static boost_matrix_float make_contingency_theorical_table(boost_matrix_float contingency_table, boost_vector_int const& _phenos_vector);
    static unsigned int sum_col(int index, boost_matrix_float contingency_table);
    static unsigned int sum_row(int index, boost_matrix_float contingency_table);
    static float compute_chi_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table);
    static unsigned int compute_liberty_degree(boost_matrix_float const& contingency_table);
    static boost_matrix_float contingency_table_conditionnal_chi_2(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector, boost_vector_string const& cond_genos_vector);
    static float compute_conditionnal_chi_2(boost_matrix_float const& contingency, unsigned int liberty_degree, vector<boost_matrix_float> & contingencies_vector);




private:

};

#endif
