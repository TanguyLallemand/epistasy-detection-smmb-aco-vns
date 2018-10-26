#include "parameters_parsing.hpp"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;


//=================================================
// Constructor
//=================================================
parameters_parsing::parameters_parsing()
{
    ifstream file("../parameters/parameters.txt");
    if(file)
    {
        string line;
        while (!file.eof())
        {
            getline(file, line);
            if (line.length() != 0 && line[0] != '#')
            {
                import_line(line);
            }
        }
    }
    else
    {
        std::cerr << "Error while opening parameters.txt !\n";
    }
}

//=================================================
// parameters_parsing : import_line
//=================================================
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
    else if(key == "alpha")
        alpha = atof(value.c_str());

    else if(key == "precision")
        precision = atof(value.c_str());

    else if(key == "subset_size_small")
        subset_size_small = atoi (value.c_str());

    else if(key == "n_trials_to_learn_mbs")
        n_trials_to_learn_mbs = atoi(value.c_str());

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

//=================================================
// parameters_parsing : split
//=================================================
vector<string> parameters_parsing::split(string const& s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim))
        tokens.push_back(item);
    return tokens;
}

//=================================================
// parameters_parsing : list_parameters
//=================================================
void parameters_parsing::list_parameters() const
{
    cout << "########### PARAMETERS ###########\n" << "header => " << header << endl
    << "separator => " << separator << endl
    << "alpha => " << alpha << endl
    << "precision => " << precision << endl
    << "aco_set_size => " << aco_set_size << endl
    << "subset_size_small => " << subset_size_small << endl
    << "n_trials_to_learn_mbs => " << n_trials_to_learn_mbs << endl
    << "n_trials_to_learn_1_mb => " << n_trials_to_learn_1_mb << endl
    << "#################################" << endl;
}

//=================================================
// parameters_parsing : update_subset_size_large
//=================================================
void parameters_parsing::update_subset_size_large(unsigned const& n_genos)
{
    if(aco_set_size == 0)
        aco_set_size = sqrt(n_genos);
}
