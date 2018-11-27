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
	int k = 3;
	std::vector<unsigned> combinations(pow(k,3));
	std::cout << combinations.size() << '\n';
	statistics::get_all_combinations(0, k, combinations);
	for (size_t i = 0; i < combinations.size(); i++) {
		std::cout << combinations[i] << '\n';
	}
    //variable_neightborhood_search vns(params);
}
