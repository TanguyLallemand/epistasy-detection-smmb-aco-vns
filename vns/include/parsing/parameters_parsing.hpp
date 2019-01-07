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
    parameters_parsing(string parameters_file);
    void import_line(std::string const& line);

// Parameters given in the OPTIONS.txt file
// Reachable from any class that include the current header (Option_file_parsing.hpp)
    int header;
    char separator;
    float alpha;

    unsigned _iteration_num;
    unsigned _pat_size_max;
    unsigned _pat_size_min;
    unsigned _max_it_vns;
    unsigned _max_it_local_search;

    string genos_file;
    string phenos_file;

    string output_directory;
    string output_prefix;



private:
    std::vector<std::string> split(std::string const& s, char delim);

};

#endif // PARAMETERS_FILE_PARSING_H
