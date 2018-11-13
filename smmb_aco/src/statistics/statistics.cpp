/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"




/*
* statistics::compute_p_value
* return p value
*/

float statistics::compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector)
{
	//TODO atention je ne redonnes pas le chi 2 score
	// Intialization of variables
	float chi_2_result = 0;
	float p_value = 0;
	unsigned int liberty_degree = 0;
	// Instanciate contigencies
	contingencies contingency_table = contingencies(2,3);
	contingencies contingency_theorical_table = contingencies(2,3);
	// Make a contingency table using datas
	contingency_table.make_contingency_table(_genos_matrix, _phenos_vector);
	contingency_theorical_table.make_contingency_theorical_table(contingency_table, _phenos_vector);
	// boost_matrix_float contingency_table = make_contingency_table(_genos_matrix, _phenos_vector);
	//TODO peut etre passer ne parametre les tailles de la contignecy table
	// Make a contingency table using datas
	// boost_matrix_float contingency_theorical_table = make_contingency_theorical_table(contingency_table, _phenos_vector);


	// Get chi square score for the two contingencies tables given
	chi_2_result = compute_chi_2(contingency_table, contingency_theorical_table);
	// Calculate liberty degree for contingency table
	liberty_degree = contingency_table.compute_liberty_degree(contingency_table);
	// Instanciate chi_squared_distribution with a given number of liberty degree
	boost::math::chi_squared_distribution<float> chi_2_distribution(liberty_degree);
	// Calculate p value following chi_squared_distribution generated and chi square score
	p_value = 1 - boost::math::cdf(chi_2_distribution, chi_2_result);
	// Return calculated p_value
	return p_value;

}



/*
* statistics::compute_chi_2
* Use a given contingency table and a theorical contingency table
* Return chi 2 score
*/

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

/*
* statistics::compute_liberty_degree
* Use a given contingency table
* Return liberty degree for this table
*/

unsigned int statistics::compute_liberty_degree(boost_matrix_float const& contingency_table)
{
	unsigned int liberty_degree = (contingency_table.size1() - 1) * (contingency_table.size2() - 1);
}


boost_matrix_float statistics::contingency_table_conditionnal_chi_2(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector, boost_vector_string const& cond_genos_vector)
{

	unsigned n_obs = _genos_matrix.size1();
	unsigned n_cond_genos = cond_genos_vector.size();
	unsigned n_contingencies = pow(3, n_cond_genos);
	unsigned int liberty_degree = 2;
	if(n_cond_genos != 0)
		liberty_degree *= 3*n_cond_genos;

	boost::numeric::ublas::vector<boost_matrix_float> contingencies_vector(n_contingencies):(2,3,0.0);
	for (size_t i = 0; i <= n_contingencies; i++) {
		contingencies_vector(i)(2,3,0.0);
	}


	//Fill contingency table (one or multiple)
	if(!cond_genos_vector.empty())
	{
		// blas::matrix_reference<blas_matrix> ref_genos_matrix = _genos_matrix.data(); // get matrix from a column
		for(unsigned i=0; i<n_obs; ++i)
		{
			// Put the current observation in the correct contingency table
			unsigned contingency_index = 0;
			unsigned j=0;
			for(list<int>::const_iterator it=cond_genos_vector.begin(); it!=cond_genos_vector.end(); ++it, ++j)
				// contingency_index += pow(3, j) * ref_genos_matrix(i, *it);
				contingency_index += pow(3, j) * _genos_matrix(i, *it);
			Contingency& c = contingencies_vector[contingency_index];
			c(_phenos_vector(i), _genos_matrix(i,0)) += 1;
		}
	}
	else
	{
		for(unsigned i=0; i<n_obs; ++i)
		{
			Contingency& c = contingencies_vector(0);
			c(_phenos_vector(i), _genos_matrix(i,0)) += 1;
		}
	}
	compute_conditionnal_chi_2(boost_matrix_float contingency, liberty_degree, contingencies_vector);
}

float statistics::compute_conditionnal_chi_2(boost_matrix_float const& contingency, unsigned int liberty_degree, vector<boost_matrix_float> & contingencies_vector)
{
	float chi_2_result = 0;
	float p_value = 0;
	for(unsigned i=0; i<contingencies_vector.size(); ++i)
	{
		G2_test_indep g2(contingencies_vector[i]);
	}
	boost::math::chi_squared_distribution<double> chi2_dist(liberty_degree);
	p_value = 1 - boost::math::cdf(chi2_dist, chi_2_result);
	if(p_value == 0)
		p_value = 2.0e-16;
}
