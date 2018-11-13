/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include "global.hpp"

/*
*
*
*
*/

    contigencies(int a, int b)
    {
        // Initialisation contingency table with floats
        boost_matrix_float _contingency_table(2,3,0.0);
    }









/*
* contigencies::make_contingency_table
* Use geno matrix and phenos vector to build a contingency table
* Return a contingency table
*/

boost_matrix_float contigencies::make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector)
{
	// For every rows
	for (size_t i = 0; i < _genos_matrix.size1(); i++)
	{
		// Store phenotype vector of current index value in a variable
		float index_row_of_contingency_table = _phenos_vector(i);
		// Store genotype matrix value of current index in a variable
		float index_col_of_contingency_table = _genos_matrix(i,0);
		// Increment contingency table at cell following index of variables given as parameters
		contingency_table.at_element(index_row_of_contingency_table, index_col_of_contingency_table) +=1;
	}
	return contingency_table;
}


/*
* contigencies::make_contingency_table
* Use geno matrix and phenos vector to build a contingency table
* Return a contingency table
*/

boost_matrix_float contigencies::make_contingency_theorical_table(boost_matrix_float contingency_table, boost_vector_int const& _phenos_vector)
{
	// // Initialisation of contingency table with floats
	// boost_matrix_float contingency_theorical_table(2,3,0.0);
	// Get vector's of phenotypes size
	unsigned int size_matrix = _phenos_vector.size();
    boost_matrix_float _contingency_theorical_table = _contingencies;


	for(unsigned i=0; i<_contingency_theorical_table.size1(); ++i)
	{
		for(unsigned j=0; j<_contingency_theorical_table.size2(); ++j)
		{
			// Theorical contingency table filling with float
			_contingency_theorical_table(i,j) = ((float)(sum_row(i,contingency_table) * (float)sum_col(j,contingency_table)) / (float)size_matrix);
		}
	}
	return _contingency_theorical_table;
}

/*
* contigencies::sum_col
* Use a given index and a contingency table
* Return sum of a given column
*/

unsigned int contigencies::sum_col(int index, boost_matrix_float const& contingency_table)
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
* contigencies::sum_row
* Use a given index and a contingency table
* Return sum of a given row
*/

unsigned int contigencies::sum_row(int index, boost_matrix_float const& contingency_table)
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
