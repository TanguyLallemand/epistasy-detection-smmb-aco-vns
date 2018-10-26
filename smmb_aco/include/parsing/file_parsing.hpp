#ifndef FILE_PARSING_HPP
#define FILE_PARSING_HPP

#include <string>
#include <boost/numeric/ublas/matrix.hpp>
#include "global.hpp"
using namespace std;

class data_parsing
{
public:
    data_parsing(string filename, int header_size, char separator);
    boost_matrix return_matrix();

private:
    string _file_name;
    int _header_size;
    int _row_number;
    int _col_number;
    char _separator;
    boost_matrix _matrix;
    void initialise_empty_matrix();
    void data_to_matrix();
    void get_col_nb(string _file_name);
    int get_line_nb(string _file_name);
};
#endif
