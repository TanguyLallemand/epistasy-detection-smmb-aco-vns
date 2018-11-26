/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "statistics.hpp"

//==============================================================================
// statistics::compute_p_value
// return p value
//==============================================================================

float statistics::compute_p_value(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector)
{
	// Intialization of variables
	float g2_result = 0;
	float p_value = 0;
	unsigned int liberty_degree = 0;

	// Container for all combinations of the current sub_subset
	list<list<unsigned>> combi_list;

	get_all_combinations(boost_vector_int & subset, list<list<unsigned>> & combi_list);
    generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset);

	// Instanciate contigencies
	contingencies contingency_table = contingencies(2,size.combi_list());
	// Make a contingency table using datas
	contingency_table.make_contingency_table(subset, _phenos_vector);
	// Make associated contingency theorical table
	contingency_table.make_contingency_theorical_table(_phenos_vector.size());
	// Get datas from contingencies class
	boost_matrix_float contingency_table_content = contingency_table.return_contingency_table();
	boost_matrix_float contingency_theorical_table_content = contingency_table.return_contingency_theorical_table();
	// Get g 2 square score for the two contingencies tables given
	g2_result = compute_g2(contingency_table_content, contingency_theorical_table_content);
	// Calculate liberty degree for contingency table
	liberty_degree = compute_liberty_degree(contingency_table_content);
	// Instanciate g_squared_distribution with a given number of liberty degree
	boost::math::chi_squared_distribution<float> g2_distribution(liberty_degree);
	// Calculate p value following g_squared_distribution generated and g square score
	p_value = 1 - boost::math::cdf(g2_distribution, g2_result);
	// Return calculated p_value
	return p_value;
}



//==============================================================================
// statistics::compute_g2
// Use a given contingency table and a theorical contingency table
// Return g 2 score
//==============================================================================

float statistics::compute_g2(boost_matrix_float const& contingency_table, boost_matrix_float const& contingency_theorical_table)
{
	float g2_result = 0;
	// Check if contingencies table are viable. In fact if one of their cell value are under 5 g 2 cannot be compute because of reliability
	if (!contingencies::reliable_test(contingency_table) || !contingencies::reliable_test(contingency_theorical_table))
	{
		return g2_result = 0.0;
	}
	// Iterate tought contingency table
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_table.size2(); ++j)
		{
			// If cell is not equal to 0
			if(contingency_table(i,j) != 0)
			{
				double div = (double) contingency_table(i,j) / contingency_theorical_table(i,j);
				g2_result += contingency_table(i,j) * log(div);
			}
		}
	}
	g2_result *= 2;
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
// smmb_aco : get_all_combinations
//==============================================================================
void smmb_aco::get_all_combinations(boost_vector_int & subset, list<list<unsigned>> & combi_list)
{
    //convert vector into list
    list<unsigned> subset(subset.begin(), subset.end());
    //temporary list that we will append to combi list every time it is modified
    list<unsigned> temp;
    //recursive function to generate all non-empty combinations
    generate_combinations(temp, combi_list, subset);
}


//==============================================================================
// smmb_aco : generate_combinations
//==============================================================================
void smmb_aco::generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset)
{
    //copy the subset for next iteration
    list<unsigned> next_subset(subset);
    //iterate the subset list
    for (auto const& h : subset)
    {
        //add current snp to the temp list
        temp.push_back(h);
        //stocking the temp in the list of combinations
        combi_list.push_back(temp);
        //remove the current x
        next_subset.remove(h);
        //recursive call on the list without current x
        generate_combinations(temp, combi_list, next_subset);
        //remove predecent snp
        temp.pop_back();
    }
}
