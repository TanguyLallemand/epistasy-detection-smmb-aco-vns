/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef CONTINGENCIES_HPP
#define CONTINGENCIES_HPP
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <vector>
#include <list>
#include "global.hpp"



class contingencies : public boost_matrix_float
{
public:
    // Constructors
    contingencies();
    contingencies(int a, int b);
    contingencies(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    contingencies(contingencies const& m);

    //setters
    void make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector);
    void make_contingency_theorical_table(int size_pheno_vector);
    // getters
    boost_matrix_float return_contingency_table();
    boost_matrix_float return_contingency_theorical_table();

    // Methods, in static to be used in conditionnal chi 2 too
    static unsigned int sum_row(int index, boost_matrix_float const& contingency_table);
    static unsigned int sum_col(int index, boost_matrix_float const& contingency_table);
    static unsigned int sum_contingency_table(boost_matrix_float const& contingency_table);
    static boost_matrix_float make_contingency_theorical_table_conditionnal(int size_pheno_vector, boost_matrix_float contingency_table);
    static std::vector<contingencies> make_contingencies_table_conditionnal(std::list<unsigned> const& cond_genos_indexes, boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost::numeric::ublas::matrix_column<boost_matrix> const& _phenos_column, int number_obs_subset, std::vector<contingencies> contingencies_vector);

    // Class Variables
    boost_matrix_float _contingency_table;
    boost_matrix_float _contingency_theorical_table;
private:

};
#endif
