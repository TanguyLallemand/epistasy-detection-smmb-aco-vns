/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef CONTINGENCIES_HPP
#define CONTINGENCIES_HPP

#include "global.hpp"


class contingencies : public boost_matrix_float
{
public:
    //==========================================================================
    // Constructors
    //==========================================================================
    // contingencies() to init a contingency table with 2,3 dimension
    contingencies();
    // contingencies(int a, int b) to init two contingencies tables with 2,3 dimension
    contingencies(int a, int b);
    // contingencies(contingencies const& m) to init a contingency table with given contingency table's dimension
    contingencies(contingencies const& m);
    //==========================================================================
    // Some static functions used in contigencies.cpp and in statistics.cpp
    //==========================================================================
    static unsigned int sum_row(int index, boost_matrix_float const& contingency_table);
    static unsigned int sum_col(int index, boost_matrix_float const& contingency_table);
    static unsigned int sum_contingency_table(boost_matrix_float const& contingency_table);
    static bool reliable_test(boost_matrix_float const& contingency_table);
    //==========================================================================
    // Some static functions used in statistics.cpp for some required operations
    // for conditionnal chi 2
    //==========================================================================
    static boost_matrix_float make_contingency_theorical_table_conditionnal( boost_matrix_float contingency_table);
    static std::vector<contingencies> make_contingencies_table_conditionnal(std::list<unsigned> const& cond_genos_indexes, boost_column const& _genos_column, boost_column const& _phenos_column, unsigned number_obs_subset, std::vector<contingencies> contingencies_vector);
    //==========================================================================
    // Class Variables
    //==========================================================================
    boost_matrix_float _contingency_table;
    boost_matrix_float _contingency_theorical_table;
};
#endif
