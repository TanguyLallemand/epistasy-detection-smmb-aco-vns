#include <iostream>
#include <stdlib.h>

#include <boost/numeric/ublas/io.hpp>

#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "smmb_aco.hpp"
#include "tools.hpp"
#include "statistics.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "global.hpp"
int main(int argc, char* argv[])
{
	// Arguments
	string genos_file = argv[1];
	string phenos_file = argv[2];

	parameters_parsing params;
	//params.list_parameters();

	params.genos_file = genos_file;
	params.phenos_file = phenos_file;
	int header = params.header;
	char separator = params.separator;

	data_parsing data(genos_file, phenos_file, header, separator);
	float pikachu;

	pikachu = statistics::compute_p_value(data._geno_matrix, data._pheno_vector);
	// Instanciation de smmb_aco
	smmb_aco test(data._geno_matrix, data._pheno_vector, params);
	test.run();

	//boost_vector mordecai = test.return_tau();
	//boost::numeric::ublas::vector<int> gg = TOOLS_HPP::sampling(params.aco_set_size, mordecai);
	//smmb_aco.run();//exemple de call de la méthode smmb. Ca ne passe pas parce qu'il faut avoir une instance d ela classe avt...


/*Algo du main

   parser les parametres
   parser les fichiers et les charger en tant que matrice
   les passer a smmb_ACO
   ecrire les resultats
 */
}
