/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>
#include <vector>

#include "statistics.hpp"
#include "contingencies.hpp"



//==============================================================================
// statistics::compute_p_value
// return p value
//==============================================================================

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
	contingency_table.make_contingency_theorical_table(_phenos_vector.size());
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
	// Return calculated p_value
	return p_value;
}



//==============================================================================
// statistics::compute_chi_2
// Use a given contingency table and a theorical contingency table
// Return chi 2 score
//==============================================================================

float statistics::compute_chi_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table)
{
	float chi_2_result = 0;
	// Check if contingencies table are viable. In fact if one of their cell value are under 5 chi 2 cannot be compute because of reliability
	if (!contingencies::reliable_test(contingency_table) || !contingencies::reliable_test(contingency_theorical_table))
	{
		return chi_2_result = 0.0;
	}
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

//==============================================================================
// statistics::compute_liberty_degree
// Use a given contingency table
// Return liberty degree for this table
//==============================================================================

unsigned int statistics::compute_liberty_degree(boost_matrix_float const& contingency_table)
{
	// Calculate liberty degree for chi 2 test, using (nbr_line-1)*(nbr_column-1)
	unsigned int liberty_degree = (contingency_table.size1() - 1) * (contingency_table.size2() - 1);
}

//==============================================================================
// statistics::make_contingencies_chi_2_conditional_test_indep
// Use a given matrix_column of genotype matrix, vector of phenotypes, and
// cond_genos_indexes
// Return chi 2 score from conditionnal chi 2
//==============================================================================

float statistics::make_contingencies_chi_2_conditional_test_indep(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost_vector_int const& _phenos_vector, std::list<unsigned> const& cond_genos_indexes)
{
	// Parse _phenos_vector to become a _phenos_column
	// Init a temporary boost matrix with size of _phenos_vector
	boost_matrix temp_pheno_matrix (_phenos_vector.size(), 1);
	for (size_t i = 0; i < _phenos_vector.size(); i++) {
		// FIll temporary matrix
		temp_pheno_matrix(i, 0) = _phenos_vector(i);
	}
	// Init a matrix column and putting temp_pheno_matrix in it
	boost::numeric::ublas::matrix_column<boost_matrix> _phenos_column (temp_pheno_matrix, 0);

	// Get some informations and stock it
	// Get number of patient in dataset
    int number_obs_subset = _genos_column.size();
	// Get number of cond_genos_indexes
    unsigned int n_cond_genos = cond_genos_indexes.size();
	// Determine how many contingencies tables are necessary
    unsigned int n_contingencies = pow(3, n_cond_genos);
	// Init liberty_degree with a default value of 2
	unsigned int liberty_degree = 2;
	// Check if liberty_degree's default is sufficient or not
    if(n_cond_genos != 0)
		// If not calculate right liberty_degree value
        liberty_degree *= 3*n_cond_genos;
	// Init a vector of contingencies table with a size of previously determined number of necessary contingency table
	std::vector<contingencies> contingencies_vector = std::vector<contingencies>(n_contingencies);
	// Build contingencies table, this function will build all contingencies table
	contingencies_vector = contingencies::make_contingencies_table_conditionnal(cond_genos_indexes, _genos_column, _phenos_column, number_obs_subset, contingencies_vector);

	//compute conditionnal chi 2, this function use contingencies table to build theorical table and will do a chi 2 test
	float result = compute_chi_2_conditional_test_indep(contingencies_vector, liberty_degree, number_obs_subset);
	// Return chi 2 scores
	return result;
}

//==============================================================================
// statistics::compute_chi_2_conditional_test_indep
// Use a vector of contigencies table, liberty degree and number of observation
// Return chi 2 score from conditionnal chi 2
//==============================================================================

float statistics::compute_chi_2_conditional_test_indep(std::vector<contingencies> contingencies_vector, unsigned int liberty_degree, unsigned int number_obs_subset)
{
	// Get number of contingencies table
	int number_contingencies = contingencies_vector.size();
	// Init p_value variable
	float p_value(number_contingencies);
	// Init chi_2_score variable
	float chi_2_score = 0;
	// Generate a chi 2 distribution for a given liberty degree
	boost::math::chi_squared_distribution<double> chi_2_distribution(liberty_degree);
	// For every contingencies tables
    for(unsigned i=0; i<number_contingencies; ++i)
    {
		// Build associated theorical table
		boost_matrix_float contingency_theorical_table_content = contingencies::make_contingency_theorical_table_conditionnal(number_obs_subset, contingencies_vector[i]);
		// Check if contingencies table are viable. In fact if one of their cell value are under 5 chi 2 cannot be compute because of reliability
		if (!contingencies::reliable_test(contingencies_vector[i]) || !contingencies::reliable_test(contingency_theorical_table_content))
		{
			chi_2_score += 0.0;
			break;
		}
		// Compute a chi 2 test
		chi_2_score += compute_chi_2(contingencies_vector[i], contingency_theorical_table_content);
    }
	// Calculate associated p value
	p_value = 1 - boost::math::cdf(chi_2_distribution, chi_2_score);
	return chi_2_score;
}
