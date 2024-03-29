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
    // contingencies(int a, int b) to init a two contingencies tables with 2,3 dimension
    contingencies(int a, int b);
    // contingencies(contingencies const& m) to init a contingency table with given contingency table's dimension
    contingencies(contingencies const& m);
    //==========================================================================
    // Setters
    //==========================================================================
    void make_contingency_table(vector<boost::numeric::ublas::matrix_column<boost_matrix>> const& pattern_datas, boost_vector_int const& _phenos_vector, vector<vector<unsigned>> all_combinations);
    void make_contingency_theorical_table(int size_pheno_vector);
    //==========================================================================
    // Getters
    //==========================================================================
    boost_matrix_float return_contingency_table();
    boost_matrix_float return_contingency_theorical_table();
    //==========================================================================
    // Some static functions used in contigencies.cpp and in statistics.cpp
    //==========================================================================
    // Compute sum of a given row
    static unsigned int sum_row(int index, boost_matrix_float const& contingency_table);
    // Compute sum of a given column
    static unsigned int sum_col(int index, boost_matrix_float const& contingency_table);
    // Compute total sm of a given contingency table
    static unsigned int sum_contingency_table(boost_matrix_float const& contingency_table);
    // Check if every cells of a given contingency table is above 5 (chi 2 application condition)
    static bool reliable_test(boost_matrix_float const& contingency_table);

    //==========================================================================
    // Class Variables
    //==========================================================================
    boost_matrix_float _contingency_table;
    boost_matrix_float _contingency_theorical_table;
};
#endif
