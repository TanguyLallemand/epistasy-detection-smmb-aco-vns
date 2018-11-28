/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "contingencies.hpp"
//==============================================================================
// Constructors
//==============================================================================
//==============================================================================
// contingencies(int a, int b) to init a two contingencies tables with
//given dimension
contingencies::contingencies(int a, int b) : boost_matrix_float(a,b)
{
    // Initialisation contingency table with floats
     _contingency_table = boost_matrix_float(a,b,0.0);
     _contingency_theorical_table = boost_matrix_float(a,b,0.0);
}
// TODO a priori a supprime seems to be useless
// // contingencies(contingencies const& m) to init a contingency table with given
// // contingency table's dimension
// contingencies::contingencies(contingencies const& m) : boost_matrix_float(m.size1(), m.size2())
// {
//     for (unsigned i = 0; i < size1(); ++i)
//     {
//         for (unsigned j = 0; j < size2(); ++j)
//         // For current contingency table, init cell(i,j) with m matrix content
//             this->at_element(i, j) = m(i,j);
//     }
// }

//==============================================================================
//Setters
//==============================================================================
//==============================================================================
//==============================================================================
// contingencies::make_contingency_table
// Use geno matrix and phenos vector to build a contingency table
// Return a contingency table
//==============================================================================

void contingencies::make_contingency_table(vector<boost::numeric::ublas::matrix_column<boost_matrix>> const& pattern_datas, boost_vector_int const& _phenos_vector, vector<vector<unsigned>> all_combinations)
{
    // For every rows
    for (size_t i = 0; i < _phenos_vector.size(); i++)
    {
        vector<unsigned> vector_pattern;
        for (size_t f = 0; f < pattern_datas.size(); f++)
        {
            vector_pattern.push_back(pattern_datas[i](f));
        }

        // Searching for vector_pattern in all_combinations
        auto it = find(all_combinations.begin(), all_combinations.end(), vector_pattern);
        // If find return a result
        if(all_combinations.end() != it)
        {
            auto index = std::distance(all_combinations.begin(), it);
            // Increment contingency table at cell following index of variables given as parameters: index of row is given by phenotype, index of column is given by index of found pattern in all possible combinations
            _contingency_table.at_element(_phenos_vector(i), index) +=1;
        }
    }
}

//==============================================================================
// contingencies::make_contingency_table
// Use geno matrix and phenos vector to build a contingency table
// Return a contingency table
//==============================================================================

void contingencies::make_contingency_theorical_table(int size_pheno_vector)
{
	// Iterating tought two dimensions of matrix
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
// g 2 test. Else return false
//==============================================================================

bool contingencies::reliable_test(boost_matrix_float const& contingency_table)
{
    // Iteration in contingency table
	for(unsigned i=0; i<contingency_table.size1(); ++i)
	{
	   for(unsigned j=0; j<contingency_table.size2(); ++j)
	   {
           // Check if each cell of contingency table is >5, if not return false because g 2 will not be reliable
           if(contingency_table(i,j) < 5)
			   return false;
	   }
	}
	return true;
}
