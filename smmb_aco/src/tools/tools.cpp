#include "tools.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <random>
#include <ctime>
//=================================================
// tools : sampling
//=================================================
boost::numeric::ublas::vector<int> tools::sampling(int subset_size, boost::numeric::ublas::vector<float> weight_vector) //generate a subset of SNP of _subset_size SNP according to weight_vector distribution
{
    //TODO voir si on donne la random seed en argument
    std::mt19937 rng;
    rng.seed(std::time(NULL));
    //std::default_random_engine rng; // random seed initialization
    boost::numeric::ublas::vector<int> SNP_picked (subset_size);
    int nb;
    for (size_t i = 0; i < subset_size; i++) {
        boost::random::discrete_distribution<int,float> distrib(weight_vector);
        nb = distrib(rng); //on pick 1 nbr avec la distrib obtenue
        SNP_picked (i) = nb; // on stocke le nbr dans la liste
        weight_vector (nb) = 0; //on passe le poids de celui qui est pick à 0 pour pas le repick
    }
    return SNP_picked; // return a subset of subset_size different SNP

}
