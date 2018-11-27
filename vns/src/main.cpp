#include <iostream>
#include <stdlib.h>
#include <string>

#include "parameters_parsing.hpp"
#include "vns.hpp"
#include "statistics.hpp"
#include "global.hpp"

int main(int argc, char* argv[])
{
    // Arguments
	std::string genos_file = argv[1];
	std::string phenos_file = argv[2];

    parameters_parsing params;


	unsigned size_of_pattern = 3;
	vector<unsigned> possible_values {0,1,2,3,4};
	vector<vector<unsigned> > all_combinations;

	statistics::init_combinations(size_of_pattern, all_combinations, possible_values);
	cout << "Total Combinations: " << all_combinations.size() << endl;

	for (int i=0; i < all_combinations.size(); i++)
	{
		cout << "{";
		for (int j=0; j < size_of_pattern; j++)
		{
			cout << all_combinations[i][j] << " ";
		}
		cout << "}" << endl;
	}

	return 0;
    //variable_neightborhood_search vns(params);
}
