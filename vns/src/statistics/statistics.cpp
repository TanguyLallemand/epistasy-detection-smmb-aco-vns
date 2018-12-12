/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "statistics.hpp"

// Objective of this file is to compute conditionnal g2 test allowing to determine
// if two variables, conditionnaly on a set of a variable are independant.
// This allows to determine if a given SNPs need to be included or removed
// from Markov Blanket in construction

//==============================================================================
// statistics::compute_p_value
// return an array with g2_score, p value and number of non reliable tests
//==============================================================================

 vector<float> statistics::compute_p_value(vector<boost::numeric::ublas::matrix_column<boost_matrix>> const& pattern_datas, boost_vector_int const& _phenos_vector)
{
	// Intialization of variables and initialization of vector to store results
	vector<float> results_to_return(3);
	vector<float> temp_result(2);
	unsigned int liberty_degree = 0;
	// Init structure storing all possible combinations
	vector<vector<unsigned>> all_combinations;
	// Get all combinations
	all_combinations = init_combinations(pattern_datas.size());
	// Instanciate contigencies
	contingencies contingency_table = contingencies(2,all_combinations.size());
	// Make a contingency table using datas
	contingency_table.make_contingency_table(pattern_datas, _phenos_vector, all_combinations);
	// Make associated contingency theorical table
	contingency_table.make_contingency_theorical_table(_phenos_vector.size());
	// Get datas from contingencies class
	boost_matrix_float contingency_table_content = contingency_table.return_contingency_table();
	boost_matrix_float contingency_theorical_table_content = contingency_table.return_contingency_theorical_table();
	// Get g2 score for the two contingencies tables given
	temp_result = compute_g2(contingency_table_content, contingency_theorical_table_content);
	// Parse results
    // Get g2 test score
	results_to_return[0] = temp_result[0];
    // Get number of non reliable test
	results_to_return[2] = temp_result[1];
	// Calculate liberty degree for contingency table
	liberty_degree = compute_liberty_degree(contingency_table_content);
	// Instanciate g_squared_distribution with a given number of liberty degree
	boost::math::chi_squared_distribution<float> g2_distribution(liberty_degree);
	// Calculate p value following g_squared_distribution generated and g square score
	results_to_return[1] = 1 - boost::math::cdf(g2_distribution, results_to_return[0]);
    // If p-value is too small, initialize it to 2.0e-16
	if (results_to_return[1] == 0)
	{
 		results_to_return[1] = 2.0e-16;
	}
	// Return calculated p_value
	return results_to_return;
}

//==============================================================================
// statistics::compute_g2
// Use a given contingency table and a theorical contingency table
// Return g 2 score
//==============================================================================

 vector<float> statistics::compute_g2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table)
{
    // Init vector to return
	vector<float> g2_result(2);
	// Iterate tought contingency table
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_table.size2(); ++j)
		{
			// If cell is not equal to 0
			if(contingency_table(i,j) != 0)
			{
				// Calculate g2 test score
				double div = (double) contingency_table(i,j) / contingency_theorical_table(i,j);
				g2_result[0] += contingency_table(i,j) * log(div);
				// Check if contingencies table are viable. In fact if one of their cell value are under 5 g2 test appear to be non reliable
				if (contingency_table(i,j) < 5 )
				{
					g2_result[1] += 1;
				}
			}
		}
	}
	g2_result[0] *= 2;
	return g2_result;
}

//==============================================================================
// statistics::compute_liberty_degree
// Use a given contingency table
// Return liberty degree for this table
//==============================================================================

unsigned statistics::compute_liberty_degree(boost_matrix_float const& contingency_table)
{
	// Calculate liberty degree for g 2 test, using (nbr_line-1)*(nbr_column-1)
	unsigned int liberty_degree = (contingency_table.size1() - 1) * (contingency_table.size2() - 1);
	return liberty_degree;
}

//==============================================================================
// statistics : get_all_combinations
// Generate all combination of pattern in a recursive fashion following pattern
// size
// return all combinations
//==============================================================================

vector<vector<unsigned>> statistics::init_combinations(unsigned pattern_size)
{
	// Add possible values for genotype
	vector<unsigned> _possible_values = {0, 1, 2};
	// Get size of pattern
	// initialisation of structure storing possible combinations
	vector<vector<unsigned>> _all_combinations;
	// initialisation of a temporary vector storing a genotype combination
    std::vector<unsigned> tuple;
	// Launch recursive function
	recursive_combination(pattern_size, tuple, _possible_values, _all_combinations);

	return _all_combinations;
}

//==============================================================================
// statistics : recursive_combination
// Generate a pattern, recursive fashion
//==============================================================================

void statistics::recursive_combination(unsigned _size_of_pattern, std::vector<unsigned> tuple, vector<unsigned> const& _possible_values, vector<vector<unsigned>> & _all_combinations)
{
	// Iterate tought list of possible values
    for (auto i : _possible_values)
	{
        // Push back current combination
    	tuple.push_back(i);
        // If there is still some possible combinations
		if (tuple.size()<_size_of_pattern)
        {
            //  Launch another time recursive function
			recursive_combination(_size_of_pattern, tuple, _possible_values, _all_combinations);
		}
		else
		{
            // Push all generated combination
			_all_combinations.push_back(tuple);
		}
        // Remove current combination
		tuple.pop_back();
    }
}
