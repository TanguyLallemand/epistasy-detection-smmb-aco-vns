#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"

//=================================================
// smmb_aco : constructeur
//=================================================
statistics::statistics(boost_matrix _genos_matrix, boost_vector_int _phenos_vector, parameters_parsing _params)
{
    _alpha_phero = _params.aco_alpha;
    _beta_phero = _params.aco_beta;
    _tau = boost_vector_int(_genos_matrix.size2(), _params.aco_tau_init);
    _eta = boost_vector_int(_genos_matrix.size2(), _params.aco_eta);
}


float statistics::compute_p_value()
{

}

boost_matrix statistics::make_contingency_table(boost_matrix & _genos_matrix, boost_vector_int & _phenos_vector)
{
    // Initialisation contingency table
    boost_matrix contingency_table(2,3,0);
    for (size_t i = 0; i < _genos_matrix.size1(); i++) {
        int row_of_contingency_table = _phenos_vector(i,0);
        int col_of_contingency_table = _genos_matrix(i,0);
        contingency_table.at_element(row_of_contingency_table, col_of_contingency_table) +=1;
    }
    return contingency_table;
}

boost_matrix statistics::make_contingency_theorical_table(boost_matrix & _genos_matrix, boost_vector_int & _phenos_vector)
{
    // Initialisation contingency table
    boost_matrix contingency_table(2,3,0);
    for (size_t i = 0; i < _genos_matrix.size1(); i++) {
        int row_of_contingency_table = _phenos_vector(i,0);
        int col_of_contingency_table = _genos_matrix(i,0);
        contingency_table.at_element(row_of_contingency_table, col_of_contingency_table) +=1;
    }
    return contingency_table;
}

unsigned int statistics::sum_col(int index)
{
    unsigned int sum_col = 0;

    for(it1 = begin1(); it1 != end1(); ++it1)
    {
        it2 = it1.begin() + index;
        s += *it2;
    }
    return sum_col;
}

unsigned int statistics::sum_row(int index)
{
    unsigned int sum_row = 0;
    blas_dmatrix::const_iterator1 it1 = begin1() + index;
    blas_dmatrix::const_iterator2 it2;

    for(it2 = it1.begin(); it2 != it1.end(); ++it2)
        s += *it2;

    return sum_row;
}
