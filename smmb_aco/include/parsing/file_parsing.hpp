#ifndef FILE_PARSING_HPP
#define FILE_PARSING_HPP

#include "global.hpp"

class data_parsing
{
public:
    //==========================================================================
    // Constructor
    //==========================================================================
    data_parsing(string _geno_filename, string _pheno_filename, int header_size, char separator);
    //==========================================================================
    // Method
    //==========================================================================
    boost_matrix return_matrix();
    //==========================================================================
    // Variables
    //==========================================================================
    boost_matrix _geno_matrix;
    boost_vector_int _pheno_vector;
    boost_vector_string _snp_id_vector;
    string _geno_filename;
    string _pheno_filename;
private:
    //==========================================================================
    // Variables
    //==========================================================================
    unsigned _header_size;
    char _separator;
    // Variables to store some internal informations
    int _row_number;
    int _col_number;
    //==========================================================================
    // Method
    //==========================================================================
    void parse_pheno();
    void parse_geno();
    void parse_snp_id();
    void get_col_nb();
    void get_line_nb();
};
#endif
