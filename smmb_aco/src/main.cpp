#include <iostream>
#include <stdlib.h>



#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "smmb_aco.hpp"

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
    //get_datas(string file_genotype, string file_phenotype)
    // Instanciation de smmb_aco
    //smmb_aco(genotype_matrix, phenotype_matrix, n_it, n_ants, K, n_it_n, alpha, tau_0, rau, tau, eta, alpha, beta);
    smmb_aco.run();//exemple de call de la méthode smmb. Ca ne passe pas parce qu'il faut avoir une instance d ela classe avt...


/*Algo du main

parser les parametres
parser les fichiers et les charger en tant que matrice
les passer a smmb_ACO
ecrire les resultats
*/
}
