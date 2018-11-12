/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"




//==============================================================================
// statistics::compute_p_value return p value
//==============================================================================
float statistics::compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector)
{
	// Intialization of variables
	float chi_2_result = 0;
	float p_value = 0;
	unsigned int liberty_degree = 0;
	// Make a contingency table using datas
	boost_matrix_float contingency_table = make_contingency_table(_genos_matrix, _phenos_vector);
	//TODO peut etre passer ne parametre les tailles de la contignecy table
	// Make a contingency table using datas
	boost_matrix_float contingency_theorical_table = make_contingency_theorical_table(contingency_table, _phenos_vector);

	std::cout << contingency_table << '\n';

	std::cout << contingency_theorical_table << '\n';
	chi_2_result = compute_chi_2(contingency_table, contingency_theorical_table);
	liberty_degree = compute_liberty_degree(contingency_table);
	boost::math::chi_squared_distribution<float> chi_2_distribution(liberty_degree);
	p_value = 1 - boost::math::cdf(chi_2_distribution, chi_2_result);
	return chi_2_result;

}

boost_matrix_float statistics::make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector)
{
	// Initialisation contingency table
	boost_matrix_float contingency_table(2,3,0.0);
	for (size_t i = 0; i < _genos_matrix.size1(); i++)
	{
		float row_of_contingency_table = _phenos_vector(i);
		float col_of_contingency_table = _genos_matrix(i,0);
		contingency_table.at_element(row_of_contingency_table, col_of_contingency_table) +=1;
	}
	return contingency_table;
}

boost_matrix_float statistics::make_contingency_theorical_table(boost_matrix_float contingency_table, boost_vector_int const& _phenos_vector)
{
	// Initialisation contingency table
	boost_matrix_float contingency_theorical_table(2,3,0.0);
	unsigned int size_matrix = _phenos_vector.size();

	// Expected contingency filling
	for(unsigned i=0; i<contingency_theorical_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_theorical_table.size2(); ++j)
		{
			contingency_theorical_table(i,j) = ((float)(sum_row(i,contingency_table) * (float)sum_col(j,contingency_table)) / (float)size_matrix);
		}
	}
	return contingency_theorical_table;
}

unsigned int statistics::sum_col(int index, boost_matrix_float contingency_table)
{
	unsigned int sum_col_of_contingency_table = 0;

	for (size_t i = 0; i < contingency_table.size1(); i++)
	{
		sum_col_of_contingency_table += contingency_table(i, index);
	}
	return sum_col_of_contingency_table;
}

unsigned int statistics::sum_row(int index, boost_matrix_float contingency_table)
{
	unsigned int sum_row_of_contingency_table = 0;

	for (size_t i = 0; i < contingency_table.size2(); i++)
	{
		sum_row_of_contingency_table += contingency_table(index, i);
	}
	return sum_row_of_contingency_table;
}

float statistics::compute_chi_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table)
{
	float chi_2_result = 0;
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_table.size2(); ++j)
		{
			if(contingency_table(i,j) != 0)
			{
				//TODO a verifier
				chi_2_result = (float) pow(contingency_table(i,j)-contingency_theorical_table(i,j), 2.0) / contingency_theorical_table(i,j);
				// chi_2_result += contingency_table(i,j) * log(div);
			}
		}
	}
	chi_2_result *= 2;
	return chi_2_result;
}

unsigned int statistics::compute_liberty_degree(boost_matrix_float const& contingency_table)
{
	unsigned int liberty_degree = (contingency_table.size1() - 1) * (contingency_table.size2() - 1);
}


boost_matrix_float statistics::contingency_table_conditionnal_chi_2(boost_vector_int const& geno_column, boost_vector_int const& _phenos_vector, boost_vector_int const& cond_genos_vector)//TODO quel type pour le vecteur de SNP?
{

    unsigned n_obs = genos.size();
    unsigned n_cond_genos = cond_genos_vector.size();
    unsigned n_contingencies = pow(3, n_cond_genos);
    liberty_degree = 2;
    if(n_cond_genos != 0)
        liberty_degree *= 3*n_cond_genos;

    _contingencies = vector<Contingency>(n_contingencies);

    // Fill contingency table (one or multiple)
    if(!cond_genos_indexes.empty())
    {
        blas::matrix_reference<blas_matrix> ref_genos_matrix = genos.data(); // get matrix from a column
        for(unsigned i=0; i<n_obs; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned contingency_index = 0;
            unsigned j=0;
            for(list<unsigned>::const_iterator it=cond_genos_indexes.begin(); it!=cond_genos_indexes.end(); ++it, ++j)
                contingency_index += pow(3, j) * ref_genos_matrix(i, *it);
            Contingency& c = _contingencies[contingency_index];
            c(phenos(i), genos(i)) += 1;
        }
    }
    else
    {
        for(unsigned i=0; i<n_obs; ++i)
        {
            Contingency& c = _contingencies[0];
            c(phenos, genos) += 1;
        }
    }

}

float statistics::compute_conditionnal_chi_2(boost_matrix_float const& contingency, liberty_degree)
{
    float chi_2_result = 0;
    for(unsigned i=0; i<_contingencies.size(); ++i)
    {
        G2_test_indep g2(_contingencies[i]);
    }
    boost::math::chi_squared_distribution<double> chi2_dist(liberty_degree);
    _pval = 1 - boost::math::cdf(chi2_dist, _g2);
    if(_pval == 0)
        _pval = 2.0e-16;
}
