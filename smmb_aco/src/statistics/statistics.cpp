#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"

//=================================================
// smmb_aco : constructeur //TODO work in progress
//=================================================
statistics::statistics(boost_matrix _genos_matrix, boost_matrix _phenos_matrix, parameters_parsing _params)
{
    _alpha_phero = _params.aco_alpha;
    _beta_phero = _params.aco_beta;
    _tau = boost_vector_float(_genos_matrix.size2(), _params.aco_tau_init);
    _eta = boost_vector_float(_genos_matrix.size2(), _params.aco_eta);
}


boost_matrix statistics::make_contingency_table(boost_matrix & _genos_matrix, boost_matrix & _phenos_matrix)
{
    // Initialisation contingency table
    boost_matrix contingency_table(2,3,0);
    for (size_t i = 0; i < _genos_matrix.size1(); i++) {
        int row_of_contingency_table = _phenos_matrix(i,0);
        int col_of_contingency_table = _genos_matrix(i,0);
        contingency_table.at_element(row_of_contingency_table, col_of_contingency_table) +=1;
    }
    return contingency_table;
}
