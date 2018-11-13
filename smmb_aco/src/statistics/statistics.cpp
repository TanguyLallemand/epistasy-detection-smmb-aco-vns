/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"
#include "contingencies.hpp"



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
	// Make a contingency table using datas
	contingency_table.make_contingency_table(_genos_matrix, _phenos_vector);
	// Make associated contingency theorical table
	contingency_table.make_contingency_theorical_table(_phenos_vector);
	// Get datas from contingencies class
	boost_matrix_float contingency_table_content = contingency_table.return_contingency_table();
	boost_matrix_float contingency_theorical_table_content = contingency_table.return_contingency_theorical_table();
	// Get chi square score for the two contingencies tables given
	chi_2_result = compute_chi_2(contingency_table_content, contingency_theorical_table_content);
	// Calculate liberty degree for contingency table
	liberty_degree = compute_liberty_degree(contingency_table_content);
	// Instanciate chi_squared_distribution with a given number of liberty degree
	boost::math::chi_squared_distribution<float> chi_2_distribution(liberty_degree);
	// Calculate p value following chi_squared_distribution generated and chi square score
	p_value = 1 - boost::math::cdf(chi_2_distribution, chi_2_result);
	std::cout << p_value << '\n';
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
