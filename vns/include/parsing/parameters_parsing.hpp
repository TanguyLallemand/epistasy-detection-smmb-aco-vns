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

    unsigned _n_it_max;
    unsigned _k_max;

    string genos_file;
    string phenos_file;

    string output_directory;
    string output_prefix;



private:
    std::vector<std::string> split(std::string const& s, char delim);

};

#endif // PARAMETERS_FILE_PARSING_H
