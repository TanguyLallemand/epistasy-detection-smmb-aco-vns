#ifndef VNS_HPP
#define VNS_HPP

#include "parameters_parsing.hpp"
#include "global.hpp"

class vns {

public:
    vns(parameters_parsing _params);

    void run(int x, int l_max, int k_max, int n_it_max);
private:
    int _n_it_max;
    int _k_max;

    void neighborhood_change(int x, int second_x, int k);
    void shake(int x, int k);
    void variable_neighborhood_descent(int x, int k_max);

};
#endif
