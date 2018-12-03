#include <iostream>
#include <stdlib.h>
#include <string>

#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "vns.hpp"
#include "statistics.hpp"
#include "global.hpp"

int main(int argc, char* argv[])
{
    // Arguments
	std::string genos_file = argv[1];
	std::string phenos_file = argv[2];
	std::string parameters_file = argv[3];

    parameters_parsing params(parameters_file);
	params.genos_file = genos_file;
	params.phenos_file = phenos_file;
	int header = params.header;
	char separator = params.separator;

	data_parsing data(genos_file, phenos_file, header, separator);
	vns test(data, params);
	test.run();

	return 0;
    //variable_neightborhood_search vns(params);
}
