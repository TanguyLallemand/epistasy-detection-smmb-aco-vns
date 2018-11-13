/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/

#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "global.hpp"

#ifndef CONTIGENCIES_HPP
#define CONTIGENCIES_HPP

class contingencies
{
public:
    // Constructors
    contigencies(int a, int b);
    contigencies(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);

    //
    void make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    void make_contingency_theorical_table(boost_matrix_float contingency_table, boost_vector_int const& _phenos_vector);
    // Common operations on contigencies
    unsigned int sum_row(int index, boost_matrix_float const& contingency_table);
    unsigned int sum_col(int index, boost_matrix_float const& contingency_table);



    // Variables
    boost_matrix_float _contingency_table;
    boost_matrix_float _contingency_theorical_table;
private:
};
#endif
