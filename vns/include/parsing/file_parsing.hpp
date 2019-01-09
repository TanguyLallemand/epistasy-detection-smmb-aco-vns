/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#ifndef FILE_PARSING_HPP
#define FILE_PARSING_HPP

#include "global.hpp"

class data_parsing
{
public:
    //==========================================================================
    // // Constructor
    //==========================================================================
    data_parsing(string _geno_filename, string _pheno_filename, int header_size, char separator);
    boost_matrix return_matrix();
    //==========================================================================
    // Input data
    //==========================================================================
    boost_matrix _geno_matrix;
    boost_vector_int _pheno_vector;
    boost_vector_string _snp_id_vector;

    string _geno_filename, _pheno_filename;

private:
    //==========================================================================
    // Variables
    //==========================================================================
    unsigned _header_size;
    char _separator;
    int _row_number;
    int _col_number;
    //==========================================================================
    // Methods to parse data
    //==========================================================================
    void parse_pheno();
    void parse_geno();
    void parse_snp_id();
    //==========================================================================
    // Some tools
    //==========================================================================
    void get_col_nb();
    void get_line_nb();
};
#endif
