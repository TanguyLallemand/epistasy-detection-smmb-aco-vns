#include <iostream>
#include <stdlib.h>



#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "smmb_aco.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "global.hpp"
int main(int argc, char* argv[])
{
    // Arguments
    string genos_file = argv[1];
    string phenos_file = argv[2];

    parameters_parsing params;
    params.list_parameters();

    params.genos_file = genos_file;
    params.phenos_file = phenos_file;
    int header = params.header;
    char separator = params.separator;
    data_parsing data_genos(genos_file, header, separator);
    boost_matrix matrix_genos = data_genos.return_matrix();
    data_parsing data_phenos(phenos_file, header, separator);
    boost_matrix matrix_phenos = data_phenos.return_matrix();
    // Instanciation de smmb_aco
    smmb_aco test(matrix_genos, matrix_phenos, params);
    //smmb_aco.run();//exemple de call de la méthode smmb. Ca ne passe pas parce qu'il faut avoir une instance d ela classe avt...


/*Algo du main

parser les parametres
parser les fichiers et les charger en tant que matrice
les passer a smmb_ACO
ecrire les resultats
*/
}
