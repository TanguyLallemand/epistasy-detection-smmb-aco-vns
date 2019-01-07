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
	// Arguments
	string genos_file = argv[1];
	string phenos_file = argv[2];
	string parameters_file = argv[3];

	parameters_parsing params(parameters_file);
	//params.list_parameters();

	params.genos_file = genos_file;
	params.phenos_file = phenos_file;
	int header = params.header;
	char separator = params.separator;

	data_parsing data(genos_file, phenos_file, header, separator);

	// Instanciation de smmb_aco
	smmb_aco test(data, params);
	test.run();
	return 0;
}
