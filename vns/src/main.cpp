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

    parameters_parsing params;
	params.genos_file = genos_file;
	params.phenos_file = phenos_file;
	int header = params.header;
	char separator = params.separator;

	data_parsing data(genos_file, phenos_file, header, separator);
	//maybe il prend pas le bon file parsing ça se vérifie dans le make file ça ?
	// non cest ok ca, il dis apres qu il sors du dossier vns
	//la suite aussi le prend pas en compte ;) <3
	// compilation ok, no verbose, no core dump
	// j ai une idee, parameters. txt est peut etre casse
	//sans doute, non j ai verifie...
	//je suis sur que c'est pas grand chose
	// AH ba c est sur que c la tristesse le bug... c est entre les copie colles from smmb
	vns test(data, params);
	test.run();

	return 0;
    //variable_neightborhood_search vns(params);
}
