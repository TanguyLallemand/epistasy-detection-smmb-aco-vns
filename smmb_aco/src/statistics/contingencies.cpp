/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include "global.hpp"
#include "contingencies.hpp"
/*
*
*
*
*/
contingencies::contingencies() : boost_matrix_float(2,3)
{
    for(unsigned i=0; i<size1(); ++i)
    {
        for(unsigned j=0; j<size2(); ++j)
            this->at_element(i,j) = 0;
    }
}
contingencies::contingencies(int a, int b) : boost_matrix_float(a,b)
{
    // Initialisation contingency table with floats
     _contingency_table = boost_matrix_float(a,b,0.0);
     _contingency_theorical_table = boost_matrix_float(a,b,0.0);
}

contingencies::contingencies(contingencies const& m) : boost_matrix_float(m.size1(), m.size2())
{
    for (unsigned i = 0; i < size1(); ++i)
    {
        for (unsigned j = 0; j < size2(); ++j)
            this->at_element(i, j) = m(i,j);
    }
}

/*
* contingencies::make_contingency_table
* Use geno matrix and phenos vector to build a contingency table
* Return a contingency table
*/

void contingencies::make_contingency_table(boost_matrix const& _genos_matrix, boost_vector_int const& _phenos_vector)
{
	// For every rows
	for (size_t i = 0; i < _genos_matrix.size1(); i++)
	{
		// Store phenotype vector of current index value in a variable
		float index_row_of_contingency_table = _phenos_vector(i);
		// Store genotype matrix value of current index in a variable
		float index_col_of_contingency_table = _genos_matrix(i,0);
		// Increment contingency table at cell following index of variables given as parameters
		_contingency_table.at_element(index_row_of_contingency_table, index_col_of_contingency_table) +=1;
	}
}


/*
* contingencies::make_contingency_table
* Use geno matrix and phenos vector to build a contingency table
* Return a contingency table
*/

void contingencies::make_contingency_theorical_table(int size_pheno_vector)
{
	// // Initialisation of contingency table with float
	for(unsigned i=0; i<_contingency_theorical_table.size1(); ++i)
	{
		for(unsigned j=0; j<_contingency_theorical_table.size2(); ++j)
		{
			// Theorical contingency table filling with float
			_contingency_theorical_table(i,j) = ((float)(sum_row(i,_contingency_table) * (float)sum_col(j,_contingency_table)) / (float)size_pheno_vector);
		}
	}
}

/*
* contingencies::sum_col
* Use a given index and a contingency table
* Return sum of a given column
*/

unsigned int contingencies::sum_col(int index, boost_matrix_float const& contingency_table)
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

unsigned int contingencies::sum_row(int index, boost_matrix_float const& contingency_table)
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

boost_matrix_float contingencies::return_contingency_table()
{
    return _contingency_table;
}

boost_matrix_float contingencies::return_contingency_theorical_table()
{
    return _contingency_theorical_table;
}
