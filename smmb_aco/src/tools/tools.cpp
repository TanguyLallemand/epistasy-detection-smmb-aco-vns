#include "tools.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <iostream>
#include <random>
#include <ctime>
//=================================================
// tools : sampling
//=================================================
boost::numeric::ublas::vector<int> tools::sampling(int subset_size, boost::numeric::ublas::vector<float> weight_vector, std::mt19937 & rng) //generate a subset of SNP of _subset_size SNP according to weight_vector distribution
{
	boost::numeric::ublas::vector<int> SNP_picked(subset_size, 0);
	int nb = 0;
	for (size_t i = 0; i < subset_size; i++)
	{
		std::cout << "hash" << '\n';
		boost::random::discrete_distribution<int,float> distrib(float weight_vector);
		// std::cout << distrib << '\n';
		nb = distrib(rng); //on pick 1 nbr avec la distrib obtenue
		SNP_picked (i) = nb; // on stocke le nbr dans la liste
		std::cout << SNP_picked << '\n';
		weight_vector (nb) = 0; //on passe le poids de celui qui est pick à 0 pour pas le repick
		std::cout << "/* message */" << '\n';
	}
	return SNP_picked; // return a subset of subset_size different SNP

}
