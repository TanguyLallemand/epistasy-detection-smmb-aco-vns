/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "statistics.hpp"

//==============================================================================
// statistics::compute_p_value
// return p value
//==============================================================================

 vector<float> statistics::compute_p_value(vector<boost::numeric::ublas::matrix_column<boost_matrix>> const& pattern_datas, boost_vector_int const& _phenos_vector)
{
	// Intialization of variables
	 vector<float> results(3);
	 vector<float> temp_result(2);
	float g2_result = 0;
	float p_value = 0;
	unsigned int liberty_degree = 0;

	vector<vector<unsigned>> all_combinations;
	all_combinations = init_combinations();

	// Instanciate contigencies
	contingencies contingency_table = contingencies(2,all_combinations.size());
	// Make a contingency table using datas
	contingency_table.make_contingency_table(pattern_datas, _phenos_vector, all_combinations);
	// Make associated contingency theorical table
	contingency_table.make_contingency_theorical_table(_phenos_vector.size());
	// Get datas from contingencies class
	boost_matrix_float contingency_table_content = contingency_table.return_contingency_table();
	boost_matrix_float contingency_theorical_table_content = contingency_table.return_contingency_theorical_table();
	// Get g 2 square score for the two contingencies tables given
	temp_result = compute_g2(contingency_table_content, contingency_theorical_table_content);
	// Parse results
	results[0] = temp_result[0];
	results[2] = temp_result[1];
	// Calculate liberty degree for contingency table
	liberty_degree = compute_liberty_degree(contingency_table_content);
	// Instanciate g_squared_distribution with a given number of liberty degree
	boost::math::chi_squared_distribution<float> g2_distribution(liberty_degree);
	// Calculate p value following g_squared_distribution generated and g square score
	results[1] = 1 - boost::math::cdf(g2_distribution, g2_result);
	// Return calculated p_value
	return results;
}



//==============================================================================
// statistics::compute_g2
// Use a given contingency table and a theorical contingency table
// Return g 2 score
//==============================================================================

 vector<float> statistics::compute_g2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table)
{
	 vector<float> g2_result(2);

	// Iterate tought contingency table
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_table.size2(); ++j)
		{
			// If cell is not equal to 0
			if(contingency_table(i,j) != 0)
			{
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
//==============================================================================

vector<vector<unsigned>> statistics::init_combinations()
{
	// Add possible values for genotype
	vector<unsigned> _possible_values = {0, 1, 2};
	// Get size of pattern
	unsigned _size_of_pattern = 3;
	// initialisation of structure storing possible combinations
	vector<vector<unsigned>> _all_combinations;
	// initialisation of a temporary vector storing a genotype combination
    std::vector<unsigned> tuple;
	// Launch recursive function
	recursive_combination(_size_of_pattern, tuple, _possible_values, _all_combinations);

	return _all_combinations;
}


void statistics::recursive_combination(unsigned _size_of_pattern, std::vector<unsigned> tuple, vector<unsigned> const& _possible_values, vector<vector<unsigned>> & _all_combinations)
{
	// Iterate tought list of possible values
    for (auto i : _possible_values)
	{
    	tuple.push_back(i);
		if (tuple.size()<_size_of_pattern){
			recursive_combination(_size_of_pattern, tuple, _possible_values, _all_combinations);
		}
		else
		{
			_all_combinations.push_back(tuple);
		}
		tuple.pop_back();
    }
}
