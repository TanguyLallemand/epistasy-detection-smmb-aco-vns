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


	statistics::test();

	return 0;
    //variable_neightborhood_search vns(params);
}
