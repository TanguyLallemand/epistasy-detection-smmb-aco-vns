#ifndef VNS_HPP
#define VNS_HPP

#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "global.hpp"

class vns
{

public:
    vns(data_parsing dataset, parameters_parsing _params);

    void run();
private:
    //input datas
    boost_matrix _genos_matrix;
    boost_vector_int _pheno_vector;
    boost_vector_string _snp_id;
    string _filename;

    //vector of patterns
    vector<list<unsigned>> pattern_list;
    //neighborhood map
    map<list<unsigned>, vector<vector<list<unsigned>>>> _neighborhood;

    //neighborhood initialisation
    void generate_patterns();
    void generate_patterns(list<unsigned> temp, list<unsigned> snp_list);
    void set_neighbors();


    int _n_it_max;
    int _k_max;
    int _l_max;

    void neighborhood_change(list<unsigned> x, list<unsigned> second_x, int k);
    list<unsigned> shake(list<unsigned> x);
    float variable_neighborhood_descent(list<unsigned> second_x, list<unsigned> & third_x);

};
#endif
