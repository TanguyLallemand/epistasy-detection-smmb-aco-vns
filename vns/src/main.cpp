/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
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
    // Get arguments passed to programm
	std::string genos_file = argv[1];
	std::string phenos_file = argv[2];
	std::string parameters_file = argv[3];
	// Parse parameters
    parameters_parsing params(parameters_file);
	params.genos_file = genos_file;
	params.phenos_file = phenos_file;
	int header = params.header;
	char separator = params.separator;
	// Parse data
	data_parsing data(genos_file, phenos_file, header, separator);
	// Instanciate vns class
	vns vns_instance(data, params);
	// Launch algorithm
	vns_instance.run();
	// End programm
	return 0;
}
