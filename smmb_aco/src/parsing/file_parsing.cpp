#include "file_parsing.hpp"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

//TODO voir si un destructeur est nécessaire

class data_parsing()
{
public:
    data_parsing(string filename, int header_size, char separator);
    void initialise_empty_matrix();
    void data_to_matrix();
    void get_col_nb(string _file_name);
    void get_line_nb(string _file_name);
    boost::numeric::ublas::matrix<int> return_matrix();

private:
    string _file_name;
    int _header_size;
    int _row_number;
    int _col_number;
    char _separator;
    boost::numeric::ublas::matrix<int> _matrix;
}

//=================================================
// data_parsing : //TODO constructeur ici
//=================================================
data_parsing::data_parsing(string filename, int header_size, char separator)
{
    _filename = filename;
    _header_size = header_size;
    _separator = separator;
    get_line_nb();
    get_col_nb();
    initialise_empty_matrix();
    data_to_matrix();
}
//=================================================
// data_parsing : data_to_matrix
//=================================================
void data_parsing::data_to_matrix()
{
    ifstream file(_file_name);
    string line;
    for (size_t x = 0; x < _header_size; x++) {
        getline(file, line)
    }
    for (size_t i = 0; i < _col_number; i++) {
        for (size_t j = 0; j < _row_number; j++) {
            getline(file, line, ',');
            _matrix (i, j) = line;
        }
    }
}

//=================================================
// data_parsing : initialise_empty_matrix
//=================================================
void data_parsing::initialise_empty_matrix()
{
    boost::numeric::ublas::matrix<int> data_matrix(int _line_number, int _col_number);
}

//=================================================
// data_parsing : get_line_nb
//=================================================
int data_parsing::get_line_nb(string _file_name)
{
    ifstream file(_file_name);
    string temp;
    _row_number=0;
    if (file) {
        while (getline(file, temp)) {
            _row_number++;
        }
        _row_number = _row_number - _header_size;
    }
    else
    {
        std::cout << "unable to open : "<< _file_name << '\n';
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
    int lenght = strlen(line);
    for (i = 0; i < lenght; i++) {
        if (line[i]==_separator) {
            _col_number++;
        }
    }
}

//=================================================
// data_parsing : return_matrix
//=================================================
boost::numeric::ublas::matrix<int> return_matrix()
{
    std::cout << _matrix << '\n';
    return _matrix;
}
