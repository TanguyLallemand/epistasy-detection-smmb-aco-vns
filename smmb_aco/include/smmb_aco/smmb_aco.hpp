#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP
#include "global.hpp"
#include <list>
#include "parameters_parsing.hpp"
using namespace std;

class smmb_aco
{
public:
    // constructeur
    smmb_aco(boost_matrix _genos_matrix, boost_vector_int _phenos_matrix, parameters_parsing _params);

    //fait tourner l'algo
    void run();

    //test TODO a remove
    boost_vector_float return_tau();
    //TODO voir pour un destructeur

private:
    //objets récupérés en argument
    boost_matrix _genos_matrix;
    boost_matrix _phenos_matrix;
    parameters_parsing _params;

    //variables initialisée par le constructeur à partir de params
    int _n_it; // nombre d'itérations ACO
    int _n_ant; // nb de fourmis
    int _subset_size; // taille du subset de variables echantillonnees à partir de _genos_matrix pour chaque fourmis (on echantillone donc les SNPs pas les individus)
    int _sub_subset_size; //taille d'une combinaison de variables echantillonnées parmi _subset_size
    int _n_it_n; // nombre d'itération maximales pour explorer l'espace de recherche
    float _alpha_stat; // seuil de significativité
    float _tau_0; // valeur initiale de phéromone de chaque variable, au debut des temps egalite parfaite car on a pas de connaissance. A traiter comme un vecteur
    //Mise a jour du taux de phéromones
    float _rho; // taux d'évaporation
    float _lambda; // values used in evaporation rates updates
    // Deux constantes utilisées pour ajuster les poids respectifs entre les taux de phéromones et les connaissances a priori. Y a peut etre des SNPs qu on connait et donc on lui donne une bonne note. Ca permet donc de regler le curseur entre importance des phéromones et importance des connaissances a priori
    float _alpha_phero;
    float _beta_phero;

    boost::numeric::ublas::vector<boost_vector_float> mem;
    list<unsigned> _mem_ant;
    //vecteur concernant les pheromones
    boost_vector_float _eta;
    boost_vector_float _tau;
    boost_vector_float _pheromone_distrib;

    // fonctions données par la prof
    void learn_MB(boost_vector_float & ant_subset, list<unsigned> & markov_blanket_a);
    void forward(bool & markov_blanket_modified, list<unsigned> & markov_blanket_a, boost_vector_float & ant_subset);
    void backward(bool & markov_blanket_modified ,list<unsigned> & markov_blanket_a);

    //fonctions qui pourrait rendre le code lisible et modulaire (by JON)
    void update_tau(); //add pheromone on a good SNP
    void sub_sampling(boost_vector_float & sub_subset, boost_vector_float & ant_subset); //compute sub_subset
    list<list<int>> get_all_combinations(boost_vector_int & sub_subset);
    void generate_combinations(list<int> temp, list<list<int>> combi_list, list<int> subset);
    void update_pheromon_distrib(); //update pheromons using _tau and _eta
};
#endif
