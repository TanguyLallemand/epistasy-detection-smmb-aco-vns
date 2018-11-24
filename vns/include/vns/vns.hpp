#ifndef VNS_HPP
#define VNS_HPP
#include "parameters_parsing.hpp"


class variable_neightborhood_search {

public:
    variable_neightborhood_search(parameters_parsing _params);
private:
    /* data */
    int _n_it_max;
    int _k_max;
};
#endif
