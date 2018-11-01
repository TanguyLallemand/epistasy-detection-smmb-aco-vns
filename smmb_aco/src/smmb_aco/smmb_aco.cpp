#include "smmb_aco.hpp"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/io.hpp>
#include "tools.hpp"
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



    _tau = boost_vector(_genos_matrix.size2(), _params.aco_tau_init); //normalement ça marche et ça init le vecteur au nbr de variable et à la valeur tau_0
}
//TODO remove it is juste a test
boost_vector smmb_aco::return_tau()
{
    std::cout << _tau << '\n';
    return _tau;
}
//=================================================
// smmb_aco : add_pheromon
//=================================================
void smmb_aco::add_pheromon(int SNP_pos)
{
    _tau[SNP_pos] += 1;//TODO cb on ajoute quand le truc est bon?
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

//=================================================
// smmb_aco : learn_MB
//=================================================
//Return Markov Blanket sous optimale eventuellemenet vide
list<unsigned> smmb_aco::learn_MB(list<unsigned> mem_a/*, P*/)
{
    list<unsigned> markov_blanket_a;
    bool markov_blanket_modified = true;
    int j = 0;
    forward(markov_blanket_modified, markov_blanket_a, j/*, P*/);
    return markov_blanket_a;
}

//=================================================
// smmb_aco : forward
//=================================================
void smmb_aco::forward(bool markov_blanket_modified, list<unsigned> markov_blanket_a, int j/*, P*/)
{
    while (markov_blanket_modified || (markov_blanket_a.empty() && j<_n_it_n))
    {
        markov_blanket_modified = false;
        /*
        TODO
        S = echantillone(P, _genos_matrix, k)
        s<-arg_max{score_association(s',T,_markov_blanket_a,mem_a)} //l'argument qui maximise
        */
        //if (p_valeur(s) << _alpha_stat) //TODO: Il faut une fonction pour calculer/renvoyer la p_valeur de la solution
        //{//cas de rejet de H_0
            //_markov_blanket_a = std::set_union (_markov_blanket_a, _markov_blanket_a + _markov_blanket_a.size(), S, S+S.size(), _markov_blanket_a.begin());// union de MB_a et S et retourne avec v.begin le debut du vecteur MB_a
            markov_blanket_modified = true;
            backward(markov_blanket_a);
        //}
        j++;
    }
}

//=================================================
// smmb_aco : backward
//=================================================
void smmb_aco::backward(list<unsigned> markov_blanket_a)
{
    for (size_t X = 0; X < markov_blanket_a.size(); X++) { // TODO it was markov_blanket instead of markov_blanket_a : erreur de frappe?
        //for (size_t S = 0; S < count; S++) {
        //TODO: pour toute combinaison S non_vides inclus dans MB
            //independance_test_conditionnal(X,T,S_0); //TODO: omg c est chaud ca
            float p_valeur = 0; //TODO temporaire pour voir si ça compile
            if (p_valeur>_alpha_stat) { //H_0: independance
                // MB <- MB\{x}; //a gerer en liste //TODO
                break;
            }
        //}
    }
}

//=================================================
// smmb_aco : run
//=================================================
void smmb_aco::run()
{
    // Initialization of Markov Blanket
    list<unsigned> markov_blanket_s;
    list<unsigned> markov_blanket_a;
    for (size_t i = 0; i < _n_it_n; i++)
    {
        //P <- calculer distribution, probabilité (tau, eta, _alpha_stat, beta)
        boost_matrix ant_colony (_n_ant, _subset_size);
        // For every ants
        for (size_t a = 0; a < _n_ant; a++)
        { // a parallelise
            boost::numeric::ublas::matrix_row<boost_matrix> ant (ant_colony, a);
            ant = TOOLS_HPP::sampling(_subset_size, _tau);
            // Initialization of memory
            list<unsigned> mem_a;
            // Generate Markov Blanket and stock it in a temp variable
            markov_blanket_a = learn_MB( mem_a/*, P*/);
        }
        // Initialization of final variable
        list<unsigned> mem;
        for (size_t a = 0; a < _n_ant; a++) {
            //ajouter(mem, mem_a); //TODO
            if (!markov_blanket_a.empty()) {
                markov_blanket_s.splice(markov_blanket_s.end(), markov_blanket_a); // move at the end of MB_s MB_a alternatively we can do MB_s.insert(MB_s.end(), MB_a.begin(), MB_a.end()); to copy MB_a content at MB_s end
            }
            //post traitement; //TODO
        }
        evaporate(); // TODO vérifier si c'est bon si on met ça la
    }
}
