/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
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
    // Check parameter parsed
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

    else if(key == "verbose")
        verbose = atoi(value.c_str());

    else if(key == "alpha")
        alpha = atof(value.c_str());

    else if(key == "subset_size_small")
        subset_size_small = atoi (value.c_str());

    else if(key == "n_trials_to_learn_1_mb")
        n_trials_to_learn_1_mb = atoi(value.c_str());

    else if(key == "gfile")
        genos_file = value;

    else if(key == "pfile")
        phenos_file = value;

    else if(key == "n_smmb_aco_runs")
        n_smmb_aco_runs = atoi(value.c_str());

    else if(key == "aco_n_ants")
        aco_n_ants = atoi(value.c_str());

    else if(key == "aco_set_size")
        aco_set_size = atoi(value.c_str());

    else if(key == "aco_n_iterations")
        aco_n_iterations = atoi(value.c_str());

    else if(key == "aco_tau_init")
        aco_tau_init = atof(value.c_str());

    else if(key == "aco_rho")
        aco_rho = atof(value.c_str());

    else if(key == "aco_lambda")
        aco_lambda = atof(value.c_str());

    else if(key == "aco_eta")
        aco_eta = atof(value.c_str());

    else if(key == "aco_alpha")
        aco_alpha = atof(value.c_str());

    else if(key == "aco_beta")
        aco_beta = atof(value.c_str());

    else {}

    n_mbs = aco_n_ants * aco_n_iterations;
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
