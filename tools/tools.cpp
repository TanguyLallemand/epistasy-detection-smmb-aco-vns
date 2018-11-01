#include "tools.hpp"
#include <boost/numeric/ublas/vector.hpp>
//=================================================
// smmb_aco : sampling
//=================================================
// TODO passer ça dans un autre fichier pour pouvoir l'utiliser dans vns
list<int> smmb_aco::sampling() //generate a subset of SNP of _subset_size SNP
{
    //TODO initialiser le random (peut etre à faire dans le constructeur)
    default_random_engine rng; // random seed initialization
    list<int> SNP_picked; //empty list of SNP to return IDEA maybe use an int boost vector
    boost_vector tau_copie = _tau;
    int nb;
    list<int> SNP_picked;
    for (size_t i = 0; i < _subset_size; i++) {
        boost::random::discrete_distribution<int,float> distrib(tau_copie);
        nb = distrib(rng); //on pick 1 nbr avec la distrib obtenue
        SNP_picked.push_back(nb); // on stocke le nbr dans la liste
        tau_copie (nb) = 0; //on passe le poids de celui qui est pick à 0 dans la copie pour pas le repick
    }
    return SNP_picked;

}
