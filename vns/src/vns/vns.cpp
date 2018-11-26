#include "vns.hpp"

variable_neightborhood_search::variable_neightborhood_search(parameters_parsing _params)
{
    this->_n_it_max = _params.n_it_max;
    //_k_max;
}

void variable_neightborhood_search::neighborhood_change(int x, int second_x, int k)
{

}

void variable_neightborhood_search::shake(int x, int k)
{

}

void variable_neightborhood_search::variable_neighborhood_descent(int x, int k_max) // This is the VND phase
{

}

void variable_neightborhood_search::run(int x, int l_max, int k_max, int n_it_max)
{
    int iterator = 0;
    int second_x;
    int third_x;
    while (n_it_max > iterator)
    {
        int k = 1;
        while (k != k_max)
        {
            shake(x, k);
            variable_neighborhood_descent(second_x, l_max);
            neighborhood_change(x, third_x, k);
        }
        iterator++;
    }
    //return x;
}
