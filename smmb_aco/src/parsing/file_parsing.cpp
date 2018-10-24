#include "file_parsing.hpp"
// on doit open les csv qui contiennent les données et les parser pour les utiliser dans les differentes techniques
#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
//TODO voir comment on retourne les données exactement
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

//initialise the data matrix
void initialise_matrix()
{
    line_nb = get_line_nb(file_genotype);
    col_nb = get_col_nb(file_genotype);
    matrix<int> data_matrix(int line_nb, int col_nb)

}

//get the line number to initialise the data matrix
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

//get the column number to initialise the data matrix
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
