#ifndef PARAMETERS_FILE_PARSING_H
#define PARAMETERS_FILE_PARSING_H

#include <string>
#include <vector>

class Parameters_file_parsing
{
public:
    Parameters_file_parsing();
    void import_line(std::string const& line);
    void list_parameters() const;
    void update_subset_size_large(unsigned const& n_genos);

// Parameters given in the OPTIONS.txt file
// Reachable from any class that include the current header (Option_file_parsing.hpp)
    int header;
    char separator;
    double alpha;
    double precision;

    unsigned n_smmb_aco_runs;

    unsigned aco_n_iterations;

//    unsigned subset_size_large;
    unsigned aco_n_ants;
    unsigned aco_set_size;
    unsigned subset_size_small;

    unsigned n_trials_to_learn_mbs;
    unsigned n_trials_to_learn_1_mb;

    double aco_tau_init;
    double aco_rho;
    double aco_lambda;
    double aco_eta;
    double aco_alpha;
    double aco_beta;

    unsigned n_mbs;

    std::string genos_file;
    std::string phenos_file;


private:
    std::vector<std::string> split(std::string const& s, char delim);

};

#endif // PARAMETERS_FILE_PARSING_H
