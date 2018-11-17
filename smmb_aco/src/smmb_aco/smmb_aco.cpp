#include "smmb_aco.hpp"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <map>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "tools.hpp"
#include "global.hpp"
#include "statistics.hpp"

using namespace std;

//==============================================================================
// smmb_aco : constructeur
//==============================================================================
smmb_aco::smmb_aco(boost_matrix genos_matrix, boost_vector_int pheno_vector, parameters_parsing _params)
{
    _genos_matrix = genos_matrix;
    _pheno_vector = pheno_vector;

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
    _eta = boost_vector_float(_genos_matrix.size2(), (float)_params.aco_eta);
    _tau = boost_vector_float(_genos_matrix.size2(), (float)_params.aco_tau_init);
    _pheromone_distrib = boost_vector_float(_genos_matrix.size2(), 0.0);

    update_pheromon_distrib(); //Initialization of the distribution for SNP sampling
}

//==============================================================================
// smmb_aco : update_tau
//==============================================================================
void smmb_aco::update_tau()
{
    //iterate through memory map
    for (auto it = _mem.begin(); it != _mem.end(); it++) {
        //evaporate the SNP selected
        _tau(it->first) = (1-_rho) * _tau(it->first);
        //iterate through the stat list of the SNP
        for (auto stat : it->second) {
            //add pheromon based on the score
            _tau(it->first) = _tau(it->first) + (_lambda * stat);
        }
    }

    // for (size_t i = 0; i < _tau.size(); i++) {
    //     //_tau(i) = (1-_rho) * _tau(i);
    //     for (size_t j = 0; j < _mem(i).size(); j++)
    //     {
    //         //_tau(i) = _tau(i) + (_lambda * _mem(i, j));
    //         _tau(i) = ((1-_rho) * _tau(i)) + _lambda * _mem(i)(j);
    //         //IDEA normalement la formule est bonne mais je trouve ça un peu con d'evaporer pour chaque stat dans la memoire (pour éviter ça je propose ce qui est en comment)
    //     }
    // }
    //repercuting changes on the distribution
    update_pheromon_distrib();
}

//==============================================================================
// smmb_aco : update_pheromon_distrib
//==============================================================================
void smmb_aco::update_pheromon_distrib()
{
    for (size_t i = 0; i < _pheromone_distrib.size(); i++) {
        _pheromone_distrib(i) = pow(_tau[i], _alpha_phero) + pow(_eta[i], _beta_phero);
        //QUESTION bon j'ai refait la formule mais reste à voir si le fait qu'on fasse la distrib juste en divisant par la somme est correct car dans la publi c'est chelou (peut etre je sais juste pas lire le langage math)
    }
}
//==============================================================================
// smmb_aco : learn_MB
//==============================================================================
//Return Markov Blanket sous optimale eventuellemenet vide
void smmb_aco::learn_MB(boost_vector_float & ant_subset, list<unsigned int> & markov_blanket_a)
{
    //to enter the loop on first iteration
    bool markov_blanket_modified = true;
    //counter of iteration number
    int j = 0;

    //on boucle pour générer la markov blanket
    while (markov_blanket_modified || (!(markov_blanket_a.empty()) && j<_n_it_n))
    {
        forward(markov_blanket_modified, markov_blanket_a, ant_subset);
        j++;
    }
    //backward(markov_blanket_modified, markov_blanket_a);
    //TODO voir si on en met un la finalement
}

//==============================================================================
// smmb_aco : forward //FIXME
//==============================================================================
void smmb_aco::forward(bool & markov_blanket_modified, list<unsigned int> & markov_blanket_a, boost_vector_float & ant_subset)
{
    //to break the loop if nothing modified
    markov_blanket_modified = false;

    //sub_subset container (S in the pseudocode)
    boost_vector_int sub_subset(_sub_subset_size, 0);
    //sub_sampling from ant_subset
    sub_sampling(sub_subset, ant_subset);
    //generating all the combination from the drawn sub_subset
    list<list<unsigned int>> combi_list;
    get_all_combinations(sub_subset, combi_list);

    float best_score = 0;
    list<unsigned int> best_pattern;
    //searching for the best combination based on score
    best_combination(best_pattern, combi_list, markov_blanket_a/*, _mem_ant*/);
        /*
        TODO
        s = argument qui maximise sur l'ensemble s' inclus ou égale à S (je considere toutes les combinaisons non vides possibles dans S ). Le truc qui est maximise c'est score d'association(s', _phenos_matrix, MB_fourmis, memoire_fourmis)
        */
        //if (p_valeur(s) << _alpha_stat) //TODO: Il faut une fonction pour calculer/renvoyer la p_valeur de la solution
        //{//rejet de l hypothese d'independance donc si on rejette on est en dependance ce qu on veut

            std::set_union (markov_blanket_a.begin(), markov_blanket_a.end(), sub_subset.begin(), sub_subset.end(), std::back_inserter(markov_blanket_a));// union de MB_a et S je crois que c'est bon //QUESTION Clement lui il modifie directement la blanket de la fourmis du coup je sais pas trop quoi penser de ton unified, mais bon comme je comprend pas ta ligne je touche pas pour le moment
            backward(markov_blanket_modified, markov_blanket_a);
            markov_blanket_modified = true;
        //}


}

//==============================================================================
// smmb_aco : backward
//==============================================================================
void smmb_aco::backward(bool & markov_blanket_modified, list<unsigned> & markov_blanket_a)
{
    if (markov_blanket_modified) {
        for (size_t X = 0; X < markov_blanket_a.size(); X++) {
            //TODO: pour toute combinaison S non_vides inclus dans MB
                //independance_test_conditionnal(X,T,S_0); //TODO: omg c est chaud ca
                // float p_valeur = statistics::compute_p_value(_genos_matrix, _phenos_matrix); //TODO temporaire pour voir si ça compile
                // if (p_valeur>_alpha_stat) { //H_0: independance
                //     auto i = std::find(begin(markov_blanket_a), end(markov_blanket_a), X);
                //     markov_blanket_a.erase(i);// MB <- MB\{x}; //veut dire MB prive de X en notation ensembliste
                //     break;
                // }
            //}
        }
    }
}

//==============================================================================
// smmb_aco : run
//==============================================================================
void smmb_aco::run()
{
    // Initialization of Markov Blanket
    list<unsigned> markov_blanket_s;
    boost::numeric::ublas::vector<list<unsigned>> markov_blanket_a(_n_ant);
    for (size_t i = 0; i < _n_it_n; i++)
    {
        // For every ants a parallelise : #pragma omp parallel for
        for (size_t a = 0; a < _n_ant; a++)
        {
            boost_vector_float ant_subset;
            ant_subset = tools::sampling(_subset_size, _pheromone_distrib); //This is the list of SNP sampled for this ant. and the distribution given on copy not ref
            // Initialization of memory
            std::map<std::vector<unsigned>, float> _mem_ant;
            // Generate Markov Blanket and stock it in a temp variable
            learn_MB(ant_subset, markov_blanket_a(a));
        }
        // Initialization of final variable
        std::map<std::vector<unsigned>, float> _mem;
        for (size_t a = 0; a < _n_ant; a++) {
            _mem.insert(_mem.end(), _mem_ant.begin(), _mem_ant.end()); //ajouter(_mem, mem_ant)
            if (!markov_blanket_a(a).empty())
            {
                markov_blanket_s.splice(markov_blanket_s.end(), markov_blanket_a(a)); // move at the end of MB_s MB_a alternatively we can do MB_s.insert(MB_s.end(), MB_a.begin(), MB_a.end()); to copy MB_a content at MB_s end // TODO a test
            }
            //post traitement; //TODO
        }
        update_tau();
    }
}

//==============================================================================
// smmb_aco : sub_sampling
//==============================================================================
void smmb_aco::sub_sampling(boost_vector_int & sub_subset, boost_vector_int const& ant_subset)
{
    boost_vector_int small_distrib(ant_subset.size()); //déclarer sous vecteur de proba pour les SNP de l'ant.
    //puis on recup les tau des snp du ant_subset
    for (size_t i = 0; i < ant_subset.size(); i++)
    {
        small_distrib (i) = _pheromone_distrib (ant_subset(i));
    }
    boost_vector_int temp;
    // on file le vecteur de proba a tools::sampling
    temp = tools::sampling(sub_subset.size(), small_distrib);

    //on prend les valeur de ant_subset aux indices renvoyés par tools::sampling
    for (size_t j = 0; j < temp.size(); j++)
    {
        sub_subset (j) = ant_subset(temp(j));
    }
}

//==============================================================================
// smmb_aco : get_all_combinations
//==============================================================================
void smmb_aco::get_all_combinations(boost_vector_int & sub_subset, list<list<unsigned int>> & combi_list)
{
    //convert vector into list
    list<unsigned int> subset(sub_subset.begin(), sub_subset.end());
    list<unsigned int> temp;
    generate_combinations(temp, combi_list, subset);
    // //reconvert the list of list to vector of vector
    // boost::numeric::ublas::vector<boost_vector_int> combi_vector(combi_list.size());
    // int i = 0;
    // for (list<list<int>>::iterator it=combi_list.begin(); it != combi_list.end(); ++it) {
    //     combi_vector(i).resize(it.size()); //TODO need to test this
    //     combi_vector(i).assign(it.begin(), it.end());
    //     i++;
    // }
    //TODO might need to go for a vector of vector or vector of list
}

//==============================================================================
// smmb_aco : generate_combinations
//==============================================================================
void smmb_aco::generate_combinations(list<unsigned int> temp, list<list<unsigned int>> & combi_list, list<unsigned int> subset)
{
    //copy the subset
    list<unsigned int> next_subset(subset);
    //iterate the subset list
    for (auto h : subset)
    {
        //add current snp to the temp list
        temp.push_back(h);
        //stocking the temp in the list of combinations
        combi_list.push_back(temp);
        //remove the current x
        next_subset.remove(h);
        //recursive call on the list without current x
        generate_combinations(temp, combi_list, next_subset);
        //remove predecent snp
        temp.pop_back();
    }
}

//==============================================================================
// smmb_aco : best_combination
//==============================================================================
void smmb_aco::best_combination(list<unsigned int> & best_pattern, list<list<unsigned int>> const& pattern_list, list<unsigned int> & markov_blanket_a/*, std::map<unsigned, list<float>> _mem_ant*/)
{
    //stock the current best score
    float best_score = 0;
    //iterate through the list of pattern
    for (auto current_pattern : pattern_list) {
        float score = 0;
        for (auto current_SNP : current_pattern) {
            //setting up the list of conditionnals SNPs
            list<unsigned int> conditionnal_set = current_pattern;
            markov_blanket_a.sort();
            conditionnal_set.sort();
            conditionnal_set.merge(markov_blanket_a);
            conditionnal_set.remove(current_SNP);
            //making a matrix column ref to pass to the function
            boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);
            //calculating score of the current SNP of the pattern and add it to the pattern score
            score += statistics::make_contingencies_chi_2_conditional_test_indep(mc, _pheno_vector, conditionnal_set);
        }
        if (score > best_score) {
            best_score = score;
            best_pattern = current_pattern;
        }
    }
}
