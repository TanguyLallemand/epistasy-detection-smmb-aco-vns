#include "file_parsing.hpp"
// on doit open les csv qui contiennent les données et les parser pour les utiliser dans les differentes techniques
#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
//TODO voir comment on retourne les données exactement (ou on met les identifiants etc)
//TODO voir pour le choix du séparateur dans le fichier input. pour le moment on met ',' par default
//IDEA Les séparateurs sont donnés en parametres dans le fichiers parameters, tt comme le nombre d eligne du header
void get_datas(string file_genotype, string file_phenotype)
{
    ifstream flux_genotype(file_genotype); //open genotype file provided (read only)

    if (flux_genotype) {

    }
    else
    {
        cout << "genotype file can't be opened" << file_genotype << endl;
    }

    ifstream flux_phenotype(file_phenotype);

    if (flux_phenotype) {
        get_line_nb(file_phenotype)
        get_col_nb(file_phenotype)
    }
    else
    {
        cout << "phenotype file can't be opened" << file_phenotype << endl;
    }


}

//=================================================
// get_datas : data_to_matrix
//=================================================
boost::numeric::ublas::matrix<string> data_to_matrix(empty_matrix, file_name)
{
    ifstream file(file_name);
    string line, value;
    int row_it = 0;
    int col_it = 0;
    while (getline(file, line)) {
        while (getline(line, value, ',')) {
            empty_matrix(row_it, col_it) = value;
            col_it++;
        }
        row_it++;
    }
    return empty_matrix;
}
//=================================================
// get_datas : initialise_matrix
//=================================================
void initialise_matrix()
{
    line_nb = get_line_nb(file_genotype);
    col_nb = get_col_nb(file_genotype);
    boost::numeric::ublas::matrix<string> data_matrix(string line_nb, string col_nb);
}

//=================================================
// get_datas : get_line_nb
//=================================================
int get_line_nb(string file_name)
{
    ifstream file(file_name);
    string temp;
    int line_nb=0;
    if (file) {
        while (getline(file, temp)) {
            line_nb++;
        }
    }
    return line_nb;
}

//=================================================
// get_datas : get_col_nb
//=================================================
int get_col_nb(string file_name, char separator)
{
    ifstream file(file_name);
    string line;
    int i;
    int col_nb=1;
    getline(file, line);
    int lenght = strlen(line);
    for (i = 0; i < lenght; i++) {
        if (line[i]==separator) {
            col_nb++;
        }
    }
}
