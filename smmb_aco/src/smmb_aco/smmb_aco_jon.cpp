#include "smmb_aco_jon.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/io.hpp>
#include "global.hpp"
using namespace std;

//=================================================
// smmb_aco : constructeur //TODO work in progress
//=================================================
smmb_aco::smmb_aco(boost_matrix _genos_matrix, boost_matrix _phenos_matrix, parameters_parsing _params)
{
    _n_it = _params.aco_n_iterations;
    _n_ant = _params.aco_n_ants;
    _rho = _params.aco_rho;
    _lambda = _params.aco_lambda;
    _alpha_phero = _params.aco_alpha;
    _beta_phero = _params.aco_beta;
    _alpha_stat = _params.alpha;



    _tau = boost_vector(_genos_matrix.size2(), params.aco_tau_init); //normalement ça marche et ça init le vecteur au nbr de variable et à la valeur tau_0
}

//=================================================
// smmb_aco : add_pheromon
//=================================================
void smmb_aco::add_pheromon(int SNP_pos)
{
    _tau[SNP_pos] += /*TODO cb on ajoute quand le truc est bon? */;
}

//=================================================
// smmb_aco : evaporate
//=================================================
void smmb_aco::evaporate()
{
    for (int i = 0; i < _tau.size(); i++) {
        _tau[i] -= _rho;
    }
}

//IDEA: Maybe drop in statistics module??, permit to be reusable for vns method.
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
