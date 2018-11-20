/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include <boost/numeric/ublas/io.hpp>
#include "global.hpp"
#include "contingencies.hpp"
//==============================================================================
// Constructors
//==============================================================================
//==============================================================================
// contingencies() to init a contingency table with 2,3 dimension
contingencies::contingencies() : boost_matrix_float(2,3)
{
    for(unsigned i=0; i<size1(); ++i)
    {
        for(unsigned j=0; j<size2(); ++j)
            // For current contingency table, init cell(i,j) with 0
            this->at_element(i,j) = 0;
    }
}
// contingencies(int a, int b) to init a two contingencies tables with
// 2,3 dimension
contingencies::contingencies(int a, int b) : boost_matrix_float(a,b)
{
    // Initialisation contingency table with floats
     _contingency_table = boost_matrix_float(a,b,0.0);
     _contingency_theorical_table = boost_matrix_float(a,b,0.0);
}
// contingencies(contingencies const& m) to init a contingency table with given
// contingency table's dimension
contingencies::contingencies(contingencies const& m) : boost_matrix_float(m.size1(), m.size2())
{
    for (unsigned i = 0; i < size1(); ++i)
    {
        for (unsigned j = 0; j < size2(); ++j)
        // For current contingency table, init cell(i,j) with m matrix content
            this->at_element(i, j) = m(i,j);
    }
}

//==============================================================================
//Setters
//==============================================================================
//==============================================================================
//==============================================================================
// contingencies::make_contingency_table
// Use geno matrix and phenos vector to build a contingency table
// Return a contingency table
//==============================================================================

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

//==============================================================================
// contingencies::make_contingency_table
// Use geno matrix and phenos vector to build a contingency table
// Return a contingency table
//==============================================================================

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

//==============================================================================
// Getters
//==============================================================================
//==============================================================================
//==============================================================================
// contingencies::return_contingency_table
// Return a contingency table
//==============================================================================

boost_matrix_float contingencies::return_contingency_table()
{
    return _contingency_table;
}

//==============================================================================
// contingencies::return_contingency_theorical_table
// Return a theorical contingency table
//==============================================================================

boost_matrix_float contingencies::return_contingency_theorical_table()
{
    return _contingency_theorical_table;
}

//==============================================================================
// Some static functions used in contigencies.cpp and in statistics.cpp
//==============================================================================
//==============================================================================
//==============================================================================
// contingencies::sum_col
// Use a given index and a contingency table
// Return sum of a given column
//==============================================================================

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

//==============================================================================
// contingencies::sum_row
// Use a given index and a contingency table
// Return sum of a given row
//==============================================================================

unsigned int contingencies::sum_row(int index, boost_matrix_float const& contingency_table)
{
	// Intialization of a variable to store sum
	unsigned int sum_row_of_contingency_table = 0;
	// For every column
	for (size_t i = 0; i < contingency_table.size2(); i++)
	{
		// Add cell content to a variable storing sum of row
		sum_row_of_contingency_table += contingency_table(index, i);
	}
	// Return sum of row
	return sum_row_of_contingency_table;
}

//==============================================================================
// contingencies::sum_contingency_table
// Use a given index and a contingency table
// Return sum of whole contingency table
//==============================================================================

unsigned int contingencies::sum_contingency_table(boost_matrix_float const& contingency_table)
{
	unsigned int sum_contingency_table = 0;
	for (size_t i = 0; i < contingency_table.size1(); i++)
	{
        // Add cell content to a variable storing sum of column
		sum_contingency_table += sum_row(i,contingency_table);
	}
	return sum_contingency_table;
}

//==============================================================================
// contingencies::reliable_test
// Use a contingency table
// Return true if all cells have a number above 5, witch is a requirement for
// chi 2 test. Else return false
//==============================================================================

bool contingencies::reliable_test(boost_matrix_float const& contingency_table)
{
    // Iteration in contingency table
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
	   for(unsigned j=0; j<contingency_table.size2(); ++j)
	   {
           // Check if each cell of contingency table is >5, if not return false because chi 2 will not be reliable
           if(contingency_table(i,j) < 5)
			   return false;
	   }
	}
	return true;
}

//==============================================================================
// Some static functions used in statistics.cpp for some required operations for
// conditionnal chi 2
//==============================================================================

std::vector<contingencies> contingencies::make_contingencies_table_conditionnal(std::list<unsigned> const& cond_genos_indexes, boost_column const& _genos_column, boost_column const& _phenos_column, unsigned number_obs_subset, std::vector<contingencies> contingencies_vector)
{
    // Fill contingency table (one or multiple)
	// Fill multiple contingency_table
    if(!cond_genos_indexes.empty())
    {
        boost::numeric::ublas::matrix<unsigned int> ref_genos_matrix;
		ref_genos_matrix = _genos_column.data(); // get matrix from a column
        for(unsigned i=0; i<number_obs_subset; ++i)
        {
            // Put the current observation in the correct contingency table
            // Intialization of iterator
            unsigned int contingency_index = 0;
            unsigned int j=0;
            // Iterate tought cond_genos_indexes list
            for(std::list<unsigned>::const_iterator it=cond_genos_indexes.begin(); it!=cond_genos_indexes.end(); ++it, ++j)
                // Find right contingency table
                contingency_index += pow(3, j) * ref_genos_matrix(i, *it);
                // Init contingency table at right index
				contingencies & c = contingencies_vector[contingency_index];
                // Fill contigency table
				unsigned cr = _phenos_column(i);
				unsigned cc = _genos_column(i);
				c(cr, cc) += 1;
        }
    }
	// Fill one contingency_table
    else
    {
        for(unsigned i=0; i<number_obs_subset; ++i)
        {
            // Fill contigency table
            contingencies & c = contingencies_vector[0];
            unsigned cr = _phenos_column(i);
            unsigned cc = _genos_column(i);
            c(cr, cc) += 1;
        }
    }
    return contingencies_vector;
}

//==============================================================================
// contingencies::make_contingency_theorical_table_independant
// Use a given index and a contingency table
// Return associated theorical contingency table
//==============================================================================

boost_matrix_float contingencies::make_contingency_theorical_table_conditionnal(boost_matrix_float contingency_table)
{
	// Initialisation of contingency table with float
	boost_matrix_float contingency_theorical_table(2,3,0.0);
    // Iterate tought two dimensions of contingency table
	for(unsigned i=0; i<contingency_theorical_table.size1(); ++i)
	{
		for(unsigned j=0; j<contingency_theorical_table.size2(); ++j)
		{
			// Theorical contingency table filling with float following this formule: (sum of each row*sum of each column)²/sum of total contingency table
			contingency_theorical_table(i,j) = ((float)(sum_row(i,contingency_table) * (float)sum_col(j,contingency_table)) / (float)sum_contingency_table(contingency_table));
		}
	}
	return contingency_theorical_table;
}
