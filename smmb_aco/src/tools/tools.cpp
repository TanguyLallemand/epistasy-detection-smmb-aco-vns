#include "tools.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <iostream>
#include <random>
#include <ctime>
//=================================================
// tools : sampling
//=================================================
boost_vector_int tools::sampling(int subset_size, boost_vector_float weight_vector, std::mt19937 & rng) //generate a subset of SNP of _subset_size SNP according to weight_vector distribution
{
	boost_vector_float SNP_picked(subset_size, 0);
	int nb = 0;
	for (size_t i = 0; i < subset_size; i++)
	{
		// std::cout << "hash" << '\n';
		boost::random::discrete_distribution<int,float> distrib(weight_vector);
		// std::cout << distrib << '\n';
		// std::cout << "a test hash 1" << '\n';
		nb = distrib(rng); //on pick 1 nbr avec la distrib obtenue
		// std::cout << "a test hash 2" << '\n';
		SNP_picked (i) = nb; // on stocke le nbr dans la liste
		// std::cout << "a test hash 3" << '\n';
		// std::cout << SNP_picked << '\n';
		weight_vector (nb) = 0; //on passe le poids de celui qui est pick à 0 pour pas le repick
		// std::cout << "a test hash 4" << '\n';
		// std::cout << SNP_picked (i)<< '\n';
	}
	// std::cout << "/* hash fin */" << '\n';
	return SNP_picked; // return a subset of subset_size different SNP


}
