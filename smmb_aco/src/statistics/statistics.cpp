/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"

//=================================================
// smmb_aco : constructeur
//=================================================



float statistics::compute_p_value(boost_matrix & _genos_matrix, boost_vector_int & _phenos_vector)
{
    boost_matrix_float contingency_table = make_contingency_table(_genos_matrix, _phenos_vector);
    //TODO peut etre passer ne parametre les tailles de la contignecy table
    boost_matrix_float contingency_theorical_table = make_contingency_theorical_table(contingency_table, _phenos_vector);

    std::cout << contingency_table << '\n';

    std::cout << contingency_theorical_table << '\n';
    float float_test = 0.2; //TODO a virer
    return float_test;

}

boost_matrix_float statistics::make_contingency_table(boost_matrix & _genos_matrix, boost_vector_int & _phenos_vector)
{
    // Initialisation contingency table
    boost_matrix_float contingency_table(2,3,0.0);
    for (size_t i = 0; i < _genos_matrix.size1(); i++) {
        float row_of_contingency_table = _phenos_vector(i);
        float col_of_contingency_table = _genos_matrix(i,0);
        contingency_table.at_element(row_of_contingency_table, col_of_contingency_table) +=1;
    }
    return contingency_table;
}

boost_matrix_float statistics::make_contingency_theorical_table(boost_matrix_float contingency_table, boost_vector_int & _phenos_vector)
{
    // Initialisation contingency table
    boost_matrix_float contingency_theorical_table(2,3,0.0);
    unsigned int size_matrix = _phenos_vector.size();

    // Expected contingency filling
    for(unsigned i=0; i<contingency_theorical_table.size1(); ++i)
    {
        for(unsigned j=0; j<contingency_theorical_table.size2(); ++j)
        {
            contingency_theorical_table(i,j) = ((float)(sum_row(i ,contingency_table) * (float)sum_col(j ,contingency_table)) / (float)size_matrix);
        }
    }
    return contingency_theorical_table;
}

unsigned int statistics::sum_col(int index, boost_matrix_float contingency_table)
{
    unsigned int sum_col_of_contingency_table = 0;

    for (size_t i = 0; i < contingency_table.size1(); i++) {
        sum_col_of_contingency_table += contingency_table(i, index);
    }
    return sum_col_of_contingency_table;
}

unsigned int statistics::sum_row(int index, boost_matrix_float contingency_table)
{
    unsigned int sum_row_of_contingency_table = 0;

    for (size_t i = 0; i < contingency_table.size2(); i++) {
        sum_row_of_contingency_table += contingency_table(index, i);
    }
    return sum_row_of_contingency_table;
}
