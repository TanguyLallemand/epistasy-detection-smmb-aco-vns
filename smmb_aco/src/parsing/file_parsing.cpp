#include "file_parsing.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/io.hpp>
#include "global.hpp"
using namespace std;

//TODO voir si un destructeur est nécessaire

//strlen c du C il vaut mieux utiliser str.size(), plus facile et c++ oriented
//=================================================
// data_parsing : constructeur
//=================================================
data_parsing::data_parsing(string filename, int header_size, char separator)
{
    _file_name = filename;
    _header_size = header_size;
    _separator = separator;
    get_line_nb();
    get_col_nb();
    initialise_empty_matrix();
    data_to_matrix();
    //std::cout << _matrix << '\n';
    std::cout << _matrix.size1() << '\n';
    std::cout << _matrix.size2() << '\n';
    std::cout << _matrix << '\n';
}
//=================================================
// data_parsing : data_to_matrix
//=================================================
void data_parsing::data_to_matrix()
{
    ifstream file(_file_name);
    string line;
    for (size_t x = 0; x < _header_size; x++) {
        getline(file, line);
    }
    for (size_t i = 0; i < _matrix.size1(); i++) {
        getline(file, line);
        std::istringstream ss(line);
        int r=0;
        for (size_t j = 0; j < line.size(); j++) {
            if (line[j]!=',') {
                _matrix(i, r) = line[j]-48;
                r++;
            }
        }
    }
}

//=================================================
// data_parsing : initialise_empty_matrix
//=================================================
void data_parsing::initialise_empty_matrix()
{
    _matrix.resize( _row_number, _col_number);
}

//=================================================
// data_parsing : get_line_nb
//=================================================
int data_parsing::get_line_nb()
{
    ifstream file(_file_name);
    string temp;
    getline(file, temp);
    _row_number=0;
    if (file) {
        while (getline(file, temp)) {
            _row_number++;
        }
        //_row_number = _row_number - _header_size;
    }
    else
    {
        cout << "unable to open : "<< _file_name << '\n';
    }
}

//=================================================
// data_parsing : get_col_nb
//=================================================
void data_parsing::get_col_nb()
{
    ifstream file(_file_name);
    string line;
    int i;
    _col_number=1;
    getline(file, line);
    int lenght = line.size();
    for (i = 0; i < lenght; i++) {
        if (line[i]==_separator) {
            _col_number++;
        }
    }
}

//=================================================
// data_parsing : return_matrix
//=================================================
boost_matrix data_parsing::return_matrix()
{
    boost_matrix _matrix;
    //std::cout << _matrix << '\n'; //Je crois que cout ne sais pas print ce type, il faut voir dans boost si ils ont pas un truc pour les prints
    return _matrix;
}
