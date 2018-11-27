#ifndef VNS_HPP
#define VNS_HPP

#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "global.hpp"

class vns
{

public:
    vns(data_parsing dataset, parameters_parsing _params);

    void run(int x, int l_max, int k_max, int n_it_max);
private:
    //input datas
    boost_matrix _genos_matrix;
    boost_vector_int _pheno_vector;
    boost_vector_string _snp_id;
    string _filename;

    //neighborhood map
    map<list<unsigned>, vector<list<list<unsigned>>>> _neighborhood;

    //neighborhood initialisation
    void generate_patterns();
    void generate_patterns(list<unsigned> temp, list<unsigned> snp_list);
    void set_neighbors();


    int _n_it_max;
    int _k_max;

    void neighborhood_change(int x, int second_x, int k);
    void shake(int x, int k);
    void variable_neighborhood_descent(int x, int k_max);

};
#endif
