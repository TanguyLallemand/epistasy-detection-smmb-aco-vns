/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "smmb_aco.hpp"
#include "tools.hpp"
#include "statistics.hpp"

#include "global.hpp"
int main(int argc, char* argv[])
{
	// Get arguments
	string genos_file = argv[1];
	string phenos_file = argv[2];
	string parameters_file = argv[3];
	// Parse parameters
	parameters_parsing params(parameters_file);
	params.genos_file = genos_file;
	params.phenos_file = phenos_file;
	int header = params.header;
	char separator = params.separator;
	// Parse data
	data_parsing data(genos_file, phenos_file, header, separator);
	// Instanciation of smmb_aco
	smmb_aco smmb_aco_instanciation(data, params);
	smmb_aco_instanciation.run();
	return 0;
}
