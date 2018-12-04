/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "statistics.hpp"


//==============================================================================
// statistics::compute_g_2
// Use a given contingency table and a theorical contingency table
// Return g 2 score
//==============================================================================

float statistics::compute_g_2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table)
{
	float g_2_result = 0;
	// Iterate tought contingency table
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_table.size2(); ++j)
		{
			// If cell is not equal to 0
			if(contingency_table(i,j) != 0)
			{
				double div = (double) contingency_table(i,j) / contingency_theorical_table(i,j);
				g_2_result += contingency_table(i,j) * log(div);
				std::cout << g_2_result << '\n';
			}
		}
	}
	g_2_result *= 2;
	return g_2_result;
}


//==============================================================================
// statistics::make_contingencies_g_2_conditional_test_indep
// Use a given matrix_column of genotype matrix, vector of phenotypes, and
// cond_genos_indexes
// Return g 2 score from conditionnal g 2
//==============================================================================

boost_vector_float statistics::make_contingencies_g_2_conditional_test_indep(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost_vector_int const& _phenos_vector, std::list<unsigned> const& cond_genos_indexes)
{
	// Parse _phenos_vector to become a _phenos_column
	// Init a temporary boost matrix with size of _phenos_vector
	boost_matrix temp_pheno_matrix (_phenos_vector.size(), 1);
	for (size_t i = 0; i < _phenos_vector.size(); i++) {
		// Fill temporary matrix
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
	{
		// If not calculate right liberty_degree value
        liberty_degree *= 3*n_cond_genos;
	}
	// Init a vector of contingencies table with a size of previously determined number of necessary contingency table
	std::vector<contingencies> contingencies_vector = std::vector<contingencies>(n_contingencies);
	// Build contingencies table, this function will build all contingencies table
	contingencies_vector = contingencies::make_contingencies_table_conditionnal(cond_genos_indexes, _genos_column, _phenos_column, number_obs_subset, contingencies_vector);

	//compute conditionnal g 2, this function use contingencies table to build theorical table and will do a g 2 test
	boost_vector_float result = compute_g_2_conditional_test_indep(contingencies_vector, liberty_degree);
	// Return g 2 scores
	return result;
}

//==============================================================================
// statistics::compute_g_2_conditional_test_indep
// Use a vector of contigencies table, liberty degree and number of observation
// Return g 2 score from conditionnal g 2
//==============================================================================
boost_vector_float statistics::compute_g_2_conditional_test_indep(std::vector<contingencies> contingencies_vector, unsigned int liberty_degree)
{
	// Init a vector that will store in first cell g 2 score and in second cell associated p value
	boost_vector_float results(3,0);
	// Get number of contingencies table
	unsigned number_contingencies = contingencies_vector.size();
	// Generate a g 2 distribution for a given liberty degree
	boost::math::chi_squared_distribution<double> g_2_distribution(liberty_degree);
	// For every contingencies tables
    for(unsigned i=0; i<number_contingencies; ++i)
    {
		// Build associated theorical table
		boost_matrix_float contingency_theorical_table_content = contingencies::make_contingency_theorical_table_conditionnal(contingencies_vector[i]);
		// Check if contingencies table are viable. In fact if one of their cell value are under 5 g 2 cannot be compute because of reliability
		if (!contingencies::reliable_test(contingencies_vector[i]) || !contingencies::reliable_test(contingency_theorical_table_content))
		{
			// If test is considered as not reliable count it
			results(2) += 1;
			break;
		}
		// Compute a g 2 test
		results(0) = compute_g_2(contingencies_vector[i], contingency_theorical_table_content);
    }
	// Calculate associated p value
	float p_value = 1 - boost::math::cdf(g_2_distribution, results(0));
	if (p_value == 0)
	{
 		results(1) = 2.0e-16;
	}
	else
	{
		results(1) = p_value;
	}

	return results;
}
