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



float statistics::make_contingencies_chi_2_conditional_test_indep(boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> const& _genos_column, boost_vector_int const& _phenos_vector, std::list<unsigned> const& cond_genos_indexes)
{

	boost_matrix temp_pheno_matrix (_phenos_vector.size(), 1);
	for (size_t i = 0; i < _phenos_vector.size(); i++) {
		temp_pheno_matrix(i, 0) = _phenos_vector(i);
	}
	std::cout << "jkb" << '\n';
	boost::numeric::ublas::matrix_column<boost_matrix> _phenos_column (temp_pheno_matrix, 0);

    unsigned int number_obs_subset = _genos_column.size();
    unsigned int n_cond_genos = cond_genos_indexes.size();
	//boost_vector_float p_value(n_cond_genos);
    unsigned int n_contingencies = pow(3, n_cond_genos);

	unsigned int liberty_degree = 2;
    if(n_cond_genos != 0)
        liberty_degree *= 3*n_cond_genos;

	std::vector<contingencies> contingencies_vector = std::vector<contingencies>(n_contingencies);

    // Fill contingency table (one or multiple)
    if(!cond_genos_indexes.empty())
    {
        boost::numeric::ublas::matrix<unsigned int> ref_genos_matrix;
		ref_genos_matrix = _genos_column.data(); // get matrix from a column
        for(unsigned i=0; i<number_obs_subset; ++i)
        {
            // Put the current observation in the correct contingency table
            unsigned int contingency_index = 0;
            unsigned int j=0;
            for(std::list<unsigned>::const_iterator it=cond_genos_indexes.begin(); it!=cond_genos_indexes.end(); ++it, ++j)
                contingency_index += pow(3, j) * ref_genos_matrix(i, *it);
				contingencies & c = contingencies_vector[contingency_index];
				unsigned cr = _phenos_column(i);
				unsigned cc = _genos_column(i);
				c(cr, cc) += 1;
        }
    }
    else
    {
        for(unsigned i=0; i<number_obs_subset; ++i)
        {
            contingencies & c = contingencies_vector[0];
            unsigned cr = _phenos_column(i);
            unsigned cc = _genos_column(i);
            c(cr, cc) += 1;
        }
    }

	//compute it
	float result = compute_chi_2_conditional_test_indep(contingencies_vector, liberty_degree, number_obs_subset);
	return result;
}
/*banniere*/
float statistics::compute_chi_2_conditional_test_indep(std::vector<contingencies> contingencies_vector, unsigned int liberty_degree, unsigned int number_obs_subset)
{
	int number_contingencies = contingencies_vector.size();
	float p_value(number_contingencies);
	boost::math::chi_squared_distribution<double> chi2_dist(liberty_degree);
	float chi_2_score = 0;
    for(unsigned i=0; i<contingencies_vector.size(); ++i)
    {
		std::cout << contingencies_vector[i] << '\n';
		boost_matrix_float contingency_theorical_table_content = statistics::make_contingency_theorical_table(number_obs_subset, contingencies_vector[i]);
		chi_2_score += compute_chi_2(contingencies_vector[i], contingency_theorical_table_content);
    }
	return chi_2_score;
}

boost_matrix_float statistics::make_contingency_theorical_table(int size_pheno_vector, boost_matrix_float contingency_table)
{
	boost_matrix_float contingency_theorical_table(2,3,0.0);
	// // Initialisation of contingency table with float
	for(unsigned i=0; i<contingency_theorical_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_theorical_table.size2(); ++j)
		{
			// Theorical contingency table filling with float
			contingency_theorical_table(i,j) = ((float)(statistics::sum_row(i,contingency_table) * (float)statistics::sum_col(j,contingency_table)) / (float)size_pheno_vector);
		}
	}
    std::cout << contingency_theorical_table << '\n';
	return contingency_theorical_table;
}

unsigned int statistics::sum_col(int index, boost_matrix_float const& contingency_table)
{
	// Intialization of a variable to store sum
	unsigned int sum_col_of_contingency_table = 0;
	// For every row
	for (size_t i = 0; i < contingency_table.size1(); i++)
	{
		// Add cell content to variable storing sum
		sum_col_of_contingency_table += contingency_table(i, index);
	}
	// Return sum of column
	return sum_col_of_contingency_table;
}

/*
* contingencies::sum_row
* Use a given index and a contingency table
* Return sum of a given row
*/

unsigned int statistics::sum_row(int index, boost_matrix_float const& contingency_table)
{
	// Intialization of a variable to store sum
	unsigned int sum_row_of_contingency_table = 0;
	// For every column
	for (size_t i = 0; i < contingency_table.size2(); i++)
	{
		// Add cell content to variable storing sum
		sum_row_of_contingency_table += contingency_table(index, i);
	}
	// Return sum of row
	return sum_row_of_contingency_table;
}
