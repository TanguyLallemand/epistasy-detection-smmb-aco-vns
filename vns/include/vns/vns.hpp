#ifndef VNS_HPP
#define VNS_HPP

#include "parameters_parsing.hpp"
#include "file_parsing.hpp"
#include "global.hpp"
#include "statistics.hpp"

class vns
{

public:
    vns(data_parsing dataset, parameters_parsing _params);

    void run();
private:
    //input datas
    boost_matrix _genos_matrix;
    boost_vector_int _phenos_vector;
    boost_vector_string _snp_id;
    string _filename;
    string _output_directory;
    string _output_prefix;

    //result map
    map<vector<unsigned>, vector<float>> _optimum_set;

    double _duration;
    int _n_it_max;
    int _k_max;
    int _l_max;

    std::mt19937 _rng;

    vector<unsigned> shake(vector<unsigned> pattern);
    vector<float> local_search(vector<unsigned> second_x, vector<unsigned> & third_x);

    void save_local_optimum(vector<unsigned> & x, vector<float> & x_score);
    void write_result_file();
    vector<float> test_pattern(vector<unsigned> const& pattern);
    vector<unsigned> generate_starting_pattern();

};
#endif
