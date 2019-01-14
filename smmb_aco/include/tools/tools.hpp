/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#ifndef TOOLS_HPP
#define TOOLS_HPP

#include "global.hpp"

class tools
{
public:
    // Function used to sample a set of SNP for each ant
    static boost::numeric::ublas::vector<int> sampling(int subset_size, boost::numeric::ublas::vector<float> weight_vector, std::mt19937 & rng);
};

#endif
