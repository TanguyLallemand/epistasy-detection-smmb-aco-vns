/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef CONTINGENCIES_HPP
#define CONTINGENCIES_HPP
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include "global.hpp"



class contingencies : public boost_matrix_float
{
public:
    // Constructors
    contingencies();
    contingencies(int a, int b);
    contingencies(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    contingencies(contingencies const& m);

    //
    void make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    void make_contingency_theorical_table(int size_pheno_vector);
    // Common operations on contingencies
    //TODO a optimise, maybe with this object
    boost_matrix_float return_contingency_table();
    boost_matrix_float return_contingency_theorical_table();



    // Variables
    boost_matrix_float _contingency_table;
    boost_matrix_float _contingency_theorical_table;
private:
    // Common operations on contingencies
    unsigned int sum_row(int index, boost_matrix_float const& contingency_table);
    unsigned int sum_col(int index, boost_matrix_float const& contingency_table);
};
#endif
