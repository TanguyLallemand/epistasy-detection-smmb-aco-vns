#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP

#include <list>
#include "parameters_parsing.hpp"
class smmb_aco
{
public:
    /*
    Parametres a donner, il faut donc les filer au constructeur?
    D: Matrice de genotype
    T: Vecteur de phénotype
    n_it: nombre d'itération ACO
    n_ants: nombre de fourmis
    K: Taille du sous ensemble de variable echantilloné à partir de D, pour chaque fourmis
    n_it_n: nombre d'itération maximale qui force l'exploration de l'espace de recherche dans Markov Blanket
    global_alpha: seuil pour erreur de type I global
    Parametres liés a l'ACO:
    tau_0: valeur intiale pour un tax de phéromone de chaque variable
    rau, lambda: Deux constantes utilisees pour mettre a jour lers taux de phéromone
    tau: vecteur indiquant les poids de connaissance a priori pour les variables
    eta:
    alpha, beta: Deux constantes utilisees pour ajuster les poids respectifs entre les taux de phéromones et les connaissances a priori
    */
// Constructeur
    smmb_aco(boost::numeric::ublas::matrix<int> genotype_matrix, boost::numeric::ublas::matrix<int> phenotype_matrix, int n_it, int n_ants, int K, int n_it_n, double alpha, int tau_0, int rau, int tau, int eta, int alpha, int beta);
    void learn_MB(boost::numeric::ublas::matrix<int> genotype_matrix, boost::numeric::ublas::matrix<int> phenotype_matrix, int K, int n_it_n, double alpha, int mem_a, P);
    void forward(MB_modifie, MB_a, int n_it_n, int j, P, boost::numeric::ublas::matrix<int> genotype_matrix, int K );
    void backward(MB_a, boost::numeric::ublas::matrix<int> phenotype_matrix, double alpha);
    void run();

private:

    int mem_a;
    // mettre les types     P;
    double alpha;
    // mettre les types     MB_A;
    int j;
    int K;
    bool MB_modified;
    // mettre les types     S;
    parameters_parsing _params;
};
#endif
