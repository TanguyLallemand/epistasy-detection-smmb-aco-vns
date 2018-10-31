#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP
#include "global.hpp"
class smmb_aco
{
public:
    // constructeur
    smmb_aco(boost_matrix _genotype, boost_matrix _phenotype, parameters_parsing _params);

    //fait tourner l'algo
    void run(boost_matrix _genos_matrix, boost_matrix _phenos_matrix, int _subset_size, size_t _n_it_n, size_t _n_ant, float _tau_0);

    //voir pour un destructeur

private:
    //objets récupérés en argument
    boost_matrix _genos_matrix;
    boost_matrix _phenos_matrix;
    parameters_parsing _params;

    //variables initialisée par le constructeur à partir de params
    int _n_it; // nombre d'itérations ACO
    int _n_ant; // nb de fourmis
    int _subset_size; // taille du subset créé par chaque fourmis
    int _n_it_n; // nombre d'itération maximales pour explorer l'espace de recherche
    double _alpha_stat; // seuil de significativité
    double _tau_0; // valeur de phéromone initiale
    double _alpha_phero; // ces 2 variables influencent l'ajustement du taux de pheromones
    double _beta_phero;
    double _rho; // taux d'évaporation
    double _lambda; // values used in evaporation rates updates
    boost_vector _eta; // vecteur de poids apriori a ajouter à tau, vector of weights (whose size is the nunmber of variables), to account for prior knowledge on the variables

    //variables modifiées pendant le run
    boost_vector _tau;//tau doit etre un vecteur de la taille du nombre de SNP

    // fonctions données par la prof
    list<unsigned> learn_MB(list<unsigned> mem_a/*, P*/);
    void forward(bool markov_blanket_modified, list<unsigned> markov_blanket_a, int j/*, P*/);
    void backward(list<unsigned> markov_blanket_a);

    //fonctions qui pourrait rendre le code lisible et modulaire (by JON)
    void add_pheromon(int SNP_pos); //add pheromone on a good SNP
    void evaporate(); //substract rho to all SNP pheromones
    void sampling() //pick a subset of SNP using tau as probability distribution
};
#endif
