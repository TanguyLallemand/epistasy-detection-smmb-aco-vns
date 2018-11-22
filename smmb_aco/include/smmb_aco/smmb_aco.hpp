#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP

#include "global.hpp"

#include "statistics.hpp"
#include "tools.hpp"
#include "parameters_parsing.hpp"

class smmb_aco
{
public:
    // constructeur
    smmb_aco(boost_matrix _genos_matrix, boost_vector_int _pheno_vector, parameters_parsing _params);

    //fait tourner l'algo
    void run();

    //test TODO a remove
    boost_vector_float return_tau();
    //TODO voir pour un destructeur

private:
    //objets récupérés en argument
    boost_matrix _genos_matrix;
    boost_vector_int _pheno_vector;
    parameters_parsing _params;

    //variables initialisée par le constructeur à partir de params
    unsigned _n_it; // nombre d'itérations ACO
    unsigned _n_ant; // nb de fourmis
    unsigned _subset_size; // taille du subset de variables echantillonnees à partir de _genos_matrix pour chaque fourmis (on echantillone donc les SNPs pas les individus)
    unsigned _sub_subset_size; //taille d'une combinaison de variables echantillonnées parmi _subset_size
    unsigned _n_it_n; // nombre d'itération maximales pour explorer l'espace de recherche
    float _alpha_stat; // seuil de significativité
    float _tau_0; // valeur initiale de phéromone de chaque variable, au debut des temps egalite parfaite car on a pas de connaissance. A traiter comme un vecteur
    //Mise a jour du taux de phéromones
    float _rho; // taux d'évaporation
    float _lambda; // values used in evaporation rates updates
    // Deux constantes utilisées pour ajuster les poids respectifs entre les taux de phéromones et les connaissances a priori. Y a peut etre des SNPs qu on connait et donc on lui donne une bonne note. Ca permet donc de regler le curseur entre importance des phéromones et importance des connaissances a priori
    float _alpha_phero;
    float _beta_phero;

    //rng seed
    std::mt19937 _rng;

    std::map<unsigned, list<float>> _mem;

    vector<map<unsigned, list<float>>> _mem_ant;

    boost::numeric::ublas::vector<list<unsigned>> _markov_blanket_a;

    // stocking the result of each run with pattern as key and occurence count of it
    map<vector<unsigned>, unsigned> markov_blanket_s;


    //vecteur concernant les pheromones
    boost_vector_float _eta;
    boost_vector_float _tau;
    boost_vector_float _pheromone_distrib;

    // main functions
    void learn_MB(boost_vector_float & ant_subset, list<unsigned> & markov_blanket_a, std::map<unsigned, list<float>> & mem_ant_ref);
    void forward(bool & markov_blanket_modified, list<unsigned> & markov_blanket_a, boost_vector_float & ant_subset, std::map<unsigned, list<float>> & mem_ant_ref);
    void backward(list<unsigned> & markov_blanket_a);

    //fonctions qui pourrait rendre le code lisible et modulaire (by JON)
    void update_tau(); //add pheromone on a good SNP
    void sub_sampling(boost_vector_int & sub_subset, boost_vector_int const& ant_subset, list<unsigned> markov_blanket_a); //compute sub_subset
    void get_all_combinations(boost_vector_int & sub_subset, list<list<unsigned>> & combi_list);
    void generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset);
    boost_vector_float best_combination(list<unsigned> & best_pattern, list<list<unsigned>> const& pattern_list, list<unsigned> & markov_blanket_a, std::map<unsigned, list<float>> & mem_ant_ref);
    void update_pheromon_distrib(); //update pheromons using _tau and _eta
};
#endif
