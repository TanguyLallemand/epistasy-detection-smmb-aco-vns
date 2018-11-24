#include "parameters_parsing.hpp"

//=================================================
// Constructor
//=================================================
parameters_parsing::parameters_parsing()
{
    ifstream file("smmb_aco/parameters/parameters.txt");
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

    else if(key == "gfile")
        genos_file = value;

    else if(key == "pfile")
        phenos_file = value;

    else if(key == "max_it")
        n_it_max = atof(value.c_str());

    else {}
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
// parameters_parsing : update_subset_size_large
//=================================================
void parameters_parsing::update_subset_size_large(unsigned const& n_genos)
{
    if(aco_set_size == 0)
        aco_set_size = sqrt(n_genos);
}