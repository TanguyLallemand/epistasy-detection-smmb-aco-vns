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
    _subset_size = _params.aco_set_size;
    _sub_subset_size = _params.subset_size_small;

    //vecteur concernant les pheromones
    _eta = boost_vector_float(_genos_matrix.size2(), _params.aco_eta);
    _tau = boost_vector_float(_genos_matrix.size2(), _params.aco_tau_init);
    _pheromone_distrib = boost_vector_float(_genos_matrix.size2(), 0);

    void update_pheromon_distrib(); //Initialization of the distribution for SNP sampling
}
//TODO remove it is juste a test
boost_vector_float smmb_aco::return_tau()
{
    std::cout << _tau << '\n';
    return _tau;
}
//=================================================
// smmb_aco : add_pheromon
//=================================================
void smmb_aco::add_pheromon(int SNP_pos)
{
    _tau[SNP_pos] += _lambda;
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
// smmb_aco : update_distribution
//=================================================
void smmb_aco::update_pheromon_distrib()
{
    for (size_t i = 0; i < _pheromone_distrib.size(); i++) {
        _pheromone_distrib.insert_element(i, _tau[i]*_alpha_phero + _eta[i]*_beta_phero);
        //insert_element (size_type i, const_reference t) Set element $i$ to the value t.
    }
}
//=================================================
// smmb_aco : learn_MB
//=================================================
//Return Markov Blanket sous optimale eventuellemenet vide
list<unsigned> smmb_aco::learn_MB(boost_vector_float ant_subset)
{
    list<unsigned> markov_blanket_a;
    bool markov_blanket_modified = true;
    int j = 0;
    while (markov_blanket_modified || (!(markov_blanket_a.empty()) && j<_n_it_n)) {
        markov_blanket_modified = false;
        boost_vector_float sub_subset(_sub_subset_size, 0);
        //sub_sampling from ant_subset
        sub_sampling(sub_subset, ant_subset);


        forward(markov_blanket_modified, markov_blanket_a, j, ant_subset);
        backward(markov_blanket_unified, ant_subset);
        j++;
    }
    return markov_blanket_a;
}

//=================================================
// smmb_aco : forward //FIXME
//=================================================
void smmb_aco::forward(bool markov_blanket_modified, list<unsigned> markov_blanket_a, int j, boost_vector_float ant_subset)
{
    //Initialise S to handle with subset
    boost::numeric::ublas::vector<int> S;
    list<unsigned> markov_blanket_unified;





        /*
        TODO
        s = argument qui maximise sur l'ensemble s' inclus ou égale à S (je considere toutes les combinaisons non vides possibles dans S ). Le truc qui est maximise c'est score d'association(s', _phenos_matrix, MB_fourmis, memoire_fourmis)
        */
        //if (p_valeur(s) << _alpha_stat) //TODO: Il faut une fonction pour calculer/renvoyer la p_valeur de la solution
        //{//rejet de l hypothese d'independance donc si on rejette on est en dependance ce qu on veut
            std::set_union (markov_blanket_a.begin(), markov_blanket_a.end(), S.begin(), S.end(), std::back_inserter(markov_blanket_unified));// union de MB_a et S je crois que c'est bon
            markov_blanket_modified = true;

        //}

    
}

//=================================================
// smmb_aco : backward
//=================================================
void smmb_aco::backward(list<unsigned> markov_blanket_a, boost_vector_float ant_subset)
{
    for (size_t X = 0; X < markov_blanket_a.size(); X++) {
        //TODO: pour toute combinaison S non_vides inclus dans MB
            //independance_test_conditionnal(X,T,S_0); //TODO: omg c est chaud ca
            float p_valeur = 0; //TODO temporaire pour voir si ça compile
            if (p_valeur>_alpha_stat) { //H_0: independance
                auto i = std::find(begin(markov_blanket_a), end(markov_blanket_a), X);
                markov_blanket_a.erase(i);// MB <- MB\{x}; //veut dire MB prive de X en notation ensembliste
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
    list<unsigned> markov_blanket_a;//attention jecrois que dans learn on reintialise ca...
    for (size_t i = 0; i < _n_it_n; i++)
    {
        // For every ants a parallelise : #pragma omp parallel for
        for (size_t a = 0; a < _n_ant; a++)
        {
            boost_vector_float ant_subset;
            ant_subset = TOOLS_HPP::sampling(_subset_size, _pheromone_distrib); //This is the list of SNP sampled for this ant.
            // Initialization of memory
            list<unsigned> _mem_ant;
            // Generate Markov Blanket and stock it in a temp variable
            markov_blanket_a = learn_MB(ant_subset);
        }
        // Initialization of final variable
        list<unsigned> mem;
        for (size_t a = 0; a < _n_ant; a++) {
            mem.insert(mem.end(), _mem_ant.begin(), _mem_ant.end()); //ajouter(mem, mem_ant)
            if (!markov_blanket_a.empty())
            {
                markov_blanket_s.splice(markov_blanket_s.end(), markov_blanket_a); // move at the end of MB_s MB_a alternatively we can do MB_s.insert(MB_s.end(), MB_a.begin(), MB_a.end()); to copy MB_a content at MB_s end // TODO a test
            }
            //post traitement; //TODO
        }
        evaporate();
    }
}

//=================================================
// smmb_aco : sub_sampling
//=================================================
void smmb_aco::sub_sampling(boost_vector_float & sub_subset, boost_vector_float ant_subset)
{
    // faut lui donner le sub_subset déja init il le prend par ref et le subset deja fait
    // et magie on a un sub_subset :D TODO
    boost_vector_float small_distrib(ant_subset.size()); //déclarer sous vecteur de proba pour les SNP de l'ant.
    //puis on recup les tau des snp du ant_subset
    for (size_t i = 0; i < ant_subset.size(); i++)
    {
        small_distrib (i) = _pheromone_distrib (ant_subset(i));
    }
    boost_vector_float temp;
    // on file le vecteur de proba a tools::sampling
    temp = TOOLS_HPP::sampling(sub_subset.size(), small_distrib);

    //on prend les valeur de ant_subset aux indices renvoyés par tools::sampling
    for (size_t j = 0; j < temp.size(); j++)
    {
        sub_subset (j) = ant_subset(temp(j));
    }
}
