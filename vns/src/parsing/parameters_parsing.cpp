#include "parameters_parsing.hpp"

//==============================================================================
// Constructor
//==============================================================================

parameters_parsing::parameters_parsing(string parameters_file)
{
    // Open parameter file
    ifstream file(parameters_file);
    // Init variable to stock line
    string line;
    // If file is opened
    if(file)
    {
        // While file is not finished to parsed
        while (!file.eof())
        {
            getline(file, line);
            // If line is not empty and begin with #
            if (!line.empty() && line[0] != '#')
            {
                // Determine what parameter is parsed and stock it in right variable
                import_line(line);
            }
        }
    }
    // Else throw an error
    else
    {
        std::cerr << "Error while opening "<< parameters_file << "\n";
    }
}


//==============================================================================
// parameters_parsing : import_line
//==============================================================================

void parameters_parsing::import_line(string const& line)
{
    vector<string> token = this->split(line, ' ');
    string const& key = token[0];
    string & value = token[1];

    if(key == "header")
        header = atoi(value.c_str());

    else if(key == "separator")
    {
        if(value == "\t")
            separator = '\t';
        else
            separator = value.at(0);
    }

    else if(key == "output_directory")
        output_directory = value;

    else if(key == "output_prefix")
        output_prefix = value;

    else if(key == "alpha")
        alpha = atof(value.c_str());

    else if(key == "gfile")
        genos_file = value;

    else if(key == "pfile")
        phenos_file = value;

    else if(key == "iteration_num")
        _iteration_num = atof(value.c_str());

    else if(key == "pat_size_max")
        _pat_size_max = atof(value.c_str());

    else if(key == "pat_size_min")
        _pat_size_min = atof(value.c_str());

    else if(key == "max_it_vns")
        _max_it_vns = atof(value.c_str());

    else if(key == "max_it_local_search")
        _max_it_local_search = atof(value.c_str());

    else {}
}

//==============================================================================
// parameters_parsing : split
//==============================================================================
// This functions allows to split line following a delimiter

vector<string> parameters_parsing::split(string const& s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim))
        tokens.push_back(item);
    return tokens;
}
