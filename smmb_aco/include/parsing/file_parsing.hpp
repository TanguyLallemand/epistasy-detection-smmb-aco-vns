#ifndef FILE_PARSING_HPP
#define FILE_PARSING_HPP

#include <string>
using namespace std;

class data_parsing
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
};
#endif
