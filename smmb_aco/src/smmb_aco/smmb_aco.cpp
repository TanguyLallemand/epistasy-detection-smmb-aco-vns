/*
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

MB_s <- 0
tau <- initialisé(tau_0)
pour i allant de n_it
    P <- calculer distribution, probabilité (tau, eta, alpha, beta)
    pour a allant de 1 a n_ants //a parrelisé?
        D_a<- echantillonner(P, D, KI)
        mem_a <- ensemble_vide
        MB_a <- learn_MB(D_a,T, k, N_it_n, global_alpha, mem_a, P)
    mem<-ensemble_vide
    pour a allant de 1 à n_ants
        ajouter(mem, mem_a)
        si (non_vide (MB_a))
        alors
            MB_S <- MB_S union {MB_a}
    post_traitement(MB_si)

//Return Markov Blanket sous optimale eventuellemenet vide
function learn_MB()
{
    Mb_a <- ensemble_vide
    MB_modified <- true
    j <- 0
    //Fast FORWARD
    tant que (MB_modified or empty(MB_a) and j<n_it_n)
        MB_modified <- false
        S<-echantillonné(P,D_a,k)
        s<-arg_max{score_association(s',T,MB_a,mem_a)} //l'argument qui maximise
        si(p_valeur(s)<global_alpha)//cas de rejet de H_0
        alors
            MB_a<-MB_a union s
            MB_modified <- true
            backward_phase(MB_a, T, global_alpha)
        j++
    return MB_a
}

procedure backward(MB, T, alpha)
{
    pour tout x element de MB
        pour toute combinaison S non_vides inclus dans MB
            realiser le terst d'independance entre X et T conditionnellement a S_0
            Si p_valeur>global_alpha//H_0: independance
            alors
            MB <- MB\{x} //a gerer en liste
            break
}
*/


#include <list>
#include "smmb_aco.hpp"
#include "statistics.hpp"
#include <boost/numeric/ublas/matrix.hpp>

class smmb_aco
{
    //=================================================
    // smmb_aco : learn_MB
    //=================================================
    //Return Markov Blanket sous optimale eventuellemenet vide
    void smmb_aco::learn_MB(boost::numeric::ublas::matrix<int> genotype_matrix, boost::numeric::ublas::matrix<int> phenotype_matrix, int K, int n_it_n, double alpha, int mem_a, P)
    {
        list<unsigned> MB_A = NULL;
        bool MB_modified = TRUE;
        int j = 0;
        forward();
        return MB_a;
    }
    //=================================================
    // smmb_aco : forward
    //=================================================
    void smmb_aco::forward(bool MB_modified, list<unsigned> MB_a, int n_it_n, int j, P, boost::numeric::ublas::matrix<int> genotype_matrix, int K )
    {
        while (MB_modified or MB_a == NULL and j<n_it_n)
        {
            MB_modified = false;
            /*
            S = echantillone(P, genotype_matrix, k)
            s<-arg_max{score_association(s',T,MB_a,mem_a)} //l'argument qui maximise
            */
            if (p_valeur(s) << alpha)
            {//cas de rejet de H_0
                MB_a = std::set_union (MB_A, MB_A + MB_A.size(), S, S+S.size(), v.begin());// union de MB_A et S et retourne avec v.begin le debut du vecteur MB_A
                MB_modified = TRUE;
                backward(MB_a, T, alpha);
            }
            j++;
        }
    }
    //=================================================
    // smmb_aco : backward
    //=================================================
    void smmb_aco::backward(list<unsigned> MB_a, boost::numeric::ublas::matrix<int> phenotype_matrix, double alpha)
    {
        for (size_t X = 0; X < MB.size; X++) {
            for (size_t S = 0; S < count; S++) {
                independance_test_conditionnal(X,T,S_0);
                if (p_valeur>alpha) { //H_0: independance
                    // MB <- MB\{x}; //a gerer en liste
                    break;
                }
            }
        }
    }
    //=================================================
    // smmb_aco : run
    //=================================================
    void smmb_aco::run()
    {
        // Initialization of Markov Blanket
        list<unsigned> MB_S = NULL;
        float tau = tau_0;
        for (size_t i = 0; i < n_iteration; i++)
        {
            //P <- calculer distribution, probabilité (tau, eta, alpha, beta)
            // For every ants
            for (size_t a = 0; a < n_ants; a++)
            { // a parallelise
                // genotype_matrix<- echantillonner(P, D, KI)
                // Initialization of memory
                list<unsigned> mem_a = NULL;
                // Generate Markov Blanket and stock it in a temp variable
                list<unsigned> MB_a = learn_MB(genotype_matrix,T, K, n_it_n, alpha, mem_a, P);
            }
            // Initialization of final variable
            list<unsigned> mem = NULL
            for (size_t a = 0; a < n_ants; a++) {
                //ajouter(mem, mem_a);
                if (MB_a != NULL) {
                    MB_S = MB_S + MB_a;
                }
                //post traitement;
            }
        }
    }
}
