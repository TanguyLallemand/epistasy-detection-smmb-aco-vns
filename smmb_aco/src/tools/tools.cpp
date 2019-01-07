/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "tools.hpp"

//=================================================
// tools : sampling
//=================================================
boost_vector_int tools::sampling(int subset_size, boost_vector_float weight_vector, std::mt19937 & rng) //generate a subset of SNP of _subset_size SNP according to weight_vector distribution
{
	boost_vector_float SNP_picked(subset_size, 0);
	int nb = 0;
	for (size_t i = 0; i < subset_size; i++)
	{
		boost::random::discrete_distribution<int,float> distrib(weight_vector);

		nb = distrib(rng); //on pick 1 nbr avec la distrib obtenue

		SNP_picked (i) = nb; // on stocke le nbr dans la liste

		weight_vector (nb) = 0; //on passe le poids de celui qui est pick à 0 pour pas le repick
	}
	return SNP_picked; // return a subset of subset_size different SNP
}
