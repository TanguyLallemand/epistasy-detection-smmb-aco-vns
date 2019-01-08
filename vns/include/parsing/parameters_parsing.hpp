/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#ifndef PARAMETERS_PARSING_HPP
#define PARAMETERS_PARSING_HPP

#include "global.hpp"
class parameters_parsing
{
public:
    //==========================================================================
    // Constructor
    //==========================================================================
    parameters_parsing(string parameters_file);
    //==========================================================================
    // Method
    //==========================================================================
    void import_line(std::string const& line);
    //==========================================================================
    // Parameters variable reachable by every class
    //==========================================================================
    // Parsing instruction
    int header;
    char separator;
    // Output parameters
    string output_directory;
    string output_prefix;
    
    float alpha;
    unsigned _iteration_num;
    unsigned _pat_size_max;
    unsigned _pat_size_min;
    unsigned _max_it_vns;
    unsigned _max_it_local_search;
    // Files given as arguments
    string genos_file;
    string phenos_file;

private:
    //==========================================================================
    // Method
    //==========================================================================
    std::vector<std::string> split(std::string const& s, char delim);

};

#endif // PARAMETERS_FILE_PARSING_H
