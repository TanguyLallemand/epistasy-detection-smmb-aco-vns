#ifndef FILE_PARSING_H
#define FILE_PARSING_H

#include <string>
using namespace std;

class file_parsing
{
    void get_datas(string file_genotype, string file_phenotype);
    void initialise_matrix();
    int get_line_nb(string file_name);
    int get_col_nb(string file_name, char separator);
}
