#include "file_parsing.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
//TODO voir comment on retourne les données exactement (ou on met les identifiants etc)
//TODO faire le constructeur
//TODO voir pour le type de la matrice car je comprends pas comment il fait le bon clement donc pour le moment j'ai mis du string pour que les identifiant passent, on verra bien à la compil ^^ #prepareuranus
//TODO voir si un destructeur est nécessaire

//Euh ca c'est pas ce qu il faut mettre dans le hpp, les proto et les variables?
class data_parsing()
{
public:
    data_parsing(string filename, int header_size, char separator)
    void initialise_empty_matrix();
    void get_data();
    void get_col_nb(string _file_name);
    void get_line_nb(string _file_name);

private:
    string _file_name;
    int _header_size;
    int _row_number;
    int _col_number;
    char _separator;
    boost::numeric::ublas::matrix<string> _matrix;
}

//=================================================
// data_parsing : //TODO constructeur ici
//=================================================

//=================================================
// data_parsing : get_data
//=================================================
void get_data()
{
    initialise_empty_matrix();
    ifstream file(_file_name);
    string line;
    for (size_t i = 0; i < _col_number; i++) {
        for (size_t j = 0; j < _row_number; j++) {
            getline(file, line, ',')
            _matrix (i, j) = line;
        }
    }
}

//=================================================
// data_parsing : initialise_empty_matrix
//=================================================
void initialise_empty_matrix()
{
    get_line_nb();
    get_col_nb();
    boost::numeric::ublas::matrix<string> data_matrix(string _line_number, string _col_number);
}

//=================================================
// data_parsing : get_line_nb
//=================================================
int get_line_nb(string _file_name)
{
    ifstream file(_file_name);
    string temp;
    _row_number=0;
    if (file) {
        while (getline(file, temp)) {
            _row_number++;
        }
    }
    else
    {
        std::cout << "unable to open : "<< _file_name << '\n';
    }
}

//=================================================
// data_parsing : get_col_nb
//=================================================
void get_col_nb(string _file_namestring _file_name)
{
    ifstream file(_file_name);
    string line;
    int i;
    _col_number=1;
    getline(file, line);
    int lenght = strlen(line);
    for (i = 0; i < lenght; i++) {
        if (line[i]==_separator) {
            _col_number++;
        }
    }
}
