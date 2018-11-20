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

    //initialisation of the markov_blanket_a to number of ant
    boost::numeric::ublas::vector<list<unsigned>> _markov_blanket_a(_n_ant);
    //initialisation of the _mem_ant to number of ant
    boost::numeric::ublas::vector<std::map<unsigned, float>> _mem_ant(_n_ant);

    update_pheromon_distrib(); //Initialization of the distribution for SNP sampling
}

//==============================================================================
// smmb_aco : update_tau
//==============================================================================
void smmb_aco::update_tau()
{
    //evaporate all the SNPs
    for (size_t i = 0; i < _tau.size(); i++) {
        _tau(i) = (1-_rho) * _tau(i);
    }
    //iterate through memory map
    for (auto it = _mem.begin(); it != _mem.end(); it++)
    {
            //add pheromon based on the score
            _tau(it->first) = _tau(it->first) + (_lambda * it->second);
    }
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
    }
}
//==============================================================================
// smmb_aco : learn_MB
//==============================================================================
//Return Markov Blanket sous optimale eventuellemenet vide
void smmb_aco::learn_MB(boost_vector_float & ant_subset, list<unsigned> & markov_blanket_a, std::map<unsigned, float> & mem_ant_ref)
{
    std::cout << "learn_MB" << '\n';
    //to enter the loop on first iteration
    bool markov_blanket_modified = true;
    //counter of iteration number
    unsigned j = 0;

    //loop to generate the markov blanket
    while (markov_blanket_modified || (!(markov_blanket_a.empty()) && j<_n_it_n))
    {
        forward(markov_blanket_modified, markov_blanket_a, ant_subset, mem_ant_ref);
        j++;
    }
    std::cout << "fin" << '\n';
    //backward(markov_blanket_modified, markov_blanket_a);
    //TODO voir si on en met un la finalement
}

//==============================================================================
// smmb_aco : forward //FIXME
//==============================================================================
void smmb_aco::forward(bool & markov_blanket_modified, list<unsigned> & markov_blanket_a, boost_vector_float & ant_subset, std::map<unsigned, float> & mem_ant_ref)
{
    std::cout << "forward" << '\n';
    //to break the loop if nothing modified
    markov_blanket_modified = false;

    // sub_subset container (S in the pseudocode)
    boost_vector_int sub_subset(_sub_subset_size, 0);
    // sub_sampling from ant_subset
    sub_sampling(sub_subset, ant_subset, markov_blanket_a);

    // generating all the combination from the drawn sub_subset
    list<list<unsigned int>> combi_list;
    get_all_combinations(sub_subset, combi_list);

    list<unsigned> best_pattern;

    // searching for the best combination based on score and returning the score and p_value of the solution
    boost_vector_float result = best_combination(best_pattern, combi_list, markov_blanket_a, mem_ant_ref);

    if (result(1) < _alpha_stat) //rejet de l hypothese d'independance donc si on rejette on est en dependance ce qu on veut
    {

        std::set_union (markov_blanket_a.begin(), markov_blanket_a.end(), best_pattern.begin(), best_pattern.end(), std::back_inserter(markov_blanket_a));// union de MB_a et S je crois que c'est bon //QUESTION Clement lui il modifie directement la blanket de la fourmis du coup je sais pas trop quoi penser de ton unified, mais bon comme je comprend pas ta ligne je touche pas pour le moment
        backward(markov_blanket_modified, markov_blanket_a);
        markov_blanket_modified = true;
    }

std::cout << "fin" << '\n';
}

//==============================================================================
// smmb_aco : backward
//==============================================================================
void smmb_aco::backward(bool & markov_blanket_modified, list<unsigned> & markov_blanket_a)
{
    std::cout << "backward" << '\n';
    list<unsigned> iterate = markov_blanket_a;
    list<unsigned> save_markov = markov_blanket_a;
    for (auto current_SNP : iterate)
    {
        list<unsigned> MB_minus_current_SNP = save_markov;
        MB_minus_current_SNP.remove(current_SNP);

        boost_vector_int vector_MB(MB_minus_current_SNP.size());
        int i=0;
        for (auto value : MB_minus_current_SNP)
        {
            vector_MB(i) = value;
            i++;
        }
        // generating all the combination from the drawn sub_subset
        list<list<unsigned>> combi_list;
        get_all_combinations(vector_MB, combi_list);
        boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);
        boost_vector_float result;
        for (auto combi : combi_list)
        {
            // Return an array with in first cell, chi 2 score and in second cell, assoicated p_value
            result = statistics::make_contingencies_chi_2_conditional_test_indep(mc, _pheno_vector, combi);
            if (result(1) > _alpha_stat)
            {
                markov_blanket_a.remove(current_SNP);
                save_markov.remove(current_SNP);
                break;
            }
        }
    }
    std::cout << "fin" << '\n';
}

//==============================================================================
// smmb_aco : run
//==============================================================================
void smmb_aco::run()
{
    // Initialization of Markov Blanket
    list<unsigned> markov_blanket_s;

    for (size_t i = 0; i < _n_it_n; i++)
    {
        //on each iteration reinitialization of ant memory and MB
        //initialisation of the markov_blanket_a to number of ant
        boost::numeric::ublas::vector<list<unsigned>> _markov_blanket_a(_n_ant);
        //initialisation of the _mem_ant to number of ant
        boost::numeric::ublas::vector<std::map<unsigned, float>> _mem_ant(_n_ant);
        //_markov_blanket_a.clear();
        //_mem_ant.clear();

        //TODO For every ants a parallelise : #pragma omp parallel for
        for (size_t a = 0; a < _n_ant; a++)
        {
            std::cout << "ant number" << '\n';
            std::cout << ' ' << a << '\n';
            boost_vector_float ant_subset;
            //This is the list of SNP sampled for this ant. and the distribution given on copy not ref to modify it
            ant_subset = tools::sampling(_subset_size, _pheromone_distrib);

            // Generate Markov Blanket
            learn_MB(ant_subset, _markov_blanket_a(a), _mem_ant(a));
        }

        _mem.clear();
        // Initialization of final variable
        for (size_t a = 0; a < _n_ant; a++)
        {
            // Insert in global map current _mem_ant
            _mem.insert(_mem_ant(a).begin(), _mem_ant(a).end());
            // If ant's markov blanket is not empty
            if (!_markov_blanket_a(a).empty())
            {
                // move at the end of MB_s MB_a alternatively we can do MB_s.insert(MB_s.end(), MB_a.begin(), MB_a.end()); to copy MB_a content at MB_s end // TODO a test
                markov_blanket_s.splice(markov_blanket_s.end(), _markov_blanket_a(a));
            }
            //post traitement; //TODO
        }

        update_tau();
    }
    std::cout << "fin" << '\n';
}

//==============================================================================
// smmb_aco : sub_sampling
//==============================================================================
void smmb_aco::sub_sampling(boost_vector_int & sub_subset, boost_vector_int const& ant_subset, list<unsigned> const& markov_blanket_a)
{
    std::cout << "sub_sampling" << '\n';
    //sub weight vector associated to the ant_subset
    boost_vector_int small_distrib(ant_subset.size());
    //getting _tau values associated to the ant_subset
    for (size_t i = 0; i < ant_subset.size(); i++)
    {
        small_distrib (i) = _pheromone_distrib (ant_subset(i));
    }
    //making markov blanket SNPs unpickable
    for (auto r : markov_blanket_a) {
        small_distrib(r) = 0;
    }

    boost_vector_int temp;
    //giving the weight vector for the ant_subset to tools::sampling
    temp = tools::sampling(sub_subset.size(), small_distrib);

    //taking the SNPs on index returned by tools::sampling in ant_subset
    for (size_t j = 0; j < temp.size(); j++)
    {
        sub_subset (j) = ant_subset(temp(j));
    }
    std::cout << "fin" << '\n';
}

//==============================================================================
// smmb_aco : get_all_combinations
//==============================================================================
void smmb_aco::get_all_combinations(boost_vector_int & sub_subset, list<list<unsigned>> & combi_list)
{
    std::cout << "get_all_combinations" << '\n';
    //convert vector into list
    list<unsigned> subset(sub_subset.begin(), sub_subset.end());
    list<unsigned> temp;
    //recursive function to generate all non-empty combinations
    generate_combinations(temp, combi_list, subset);
    std::cout << "fin" << '\n';
}

//==============================================================================
// smmb_aco : generate_combinations
//==============================================================================
void smmb_aco::generate_combinations(list<unsigned int> temp, list<list<unsigned int>> & combi_list, list<unsigned int> subset)
{
    std::cout << "generate_combinations" << '\n';
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
    std::cout << "fin" << '\n';
}

//==============================================================================
// smmb_aco : best_combination
//==============================================================================
boost_vector_float smmb_aco::best_combination(list<unsigned int> & best_pattern, list<list<unsigned int>> const& pattern_list, list<unsigned> & markov_blanket_a, std::map<unsigned, float> & mem_ant_ref)
{
    std::cout << "best_combination" << '\n';
    //stock the current best_result
    boost_vector_float best_result(2, 0);
    //iterate through the list of pattern
    for (auto current_pattern : pattern_list)
    {
        std::cout << "/* message 1*/" << '\n';
        boost_vector_float result_pattern(2, 0);
        for (auto current_SNP : current_pattern)
        {
            std::cout << "/* message 2*/" << '\n';
            //setting up the list of conditionnals SNPs
            list<unsigned int> conditionnal_set = current_pattern;
            markov_blanket_a.sort();
            conditionnal_set.sort();
            conditionnal_set.merge(markov_blanket_a);
            conditionnal_set.remove(current_SNP);
            std::cout << "/* message 3*/" << '\n';
            //making a matrix column ref to pass to the function
            std::cout << current_SNP << '\n';
            boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);

            std::cout << mc << '\n';
            std::cout << "/* message 3.5*/" << '\n';
            //calculating score of the current SNP of the pattern and add it to the pattern score
            boost_vector_float result_SNP = statistics::make_contingencies_chi_2_conditional_test_indep(mc, _pheno_vector, conditionnal_set);
            std::cout << "/* message 4*/" << '\n';
            //stocking result in the ant_memory
            mem_ant_ref.insert(std::pair<unsigned, float>(current_SNP, result_SNP(0)));

            //score of the pattern
            if (result_SNP(0) > result_pattern(0))
            {
                result_pattern = result_SNP;
            }
            std::cout << "/* message 5*/" << '\n';
        }
        //if the score of the current pattern is better than the current one it become the current best pattern
        if (result_pattern(0) > best_result(0))
        {
            best_result = result_pattern;
            best_pattern = current_pattern;
        }
    }
    std::cout << "fin" << '\n';
    return best_result;
}
