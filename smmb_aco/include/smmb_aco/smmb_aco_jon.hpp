#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP
#include "global.hpp"
class smmb_aco
{
public:
    // constructeur
    smmb_aco(boost_matrix genotype, boost_matrix phenotype, parameters_parsing params);

    //fait tourner l'algo
    void run();

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
    double _lambda; // jsp
    boost_vector eta; // vecteur de poids apriori a ajouter à tau

    //variables modifiées pendant le run
    boost_vector _tau;//tau doit etre un vecteur de la taille du nombre de SNP

    // fonctions données par la prof
    list<unsigned> learn_MB(); // TODO mettre les parametres
    void forward(); // TODO mettre les parametres
    void backward(); // TODO mettre les parametres

    //fonctions qui pourrait rendre le code lisible et modulaire (by JON)
    void add_pheromon(int SNP_pos);
    void evaporate();
};
#endif
