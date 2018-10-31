_n_ant/*
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
#include "global.hpp"


//=================================================
// smmb_aco : learn_MB
//=================================================
//Return Markov Blanket sous optimale eventuellemenet vide
list<unsigned> smmb_aco::learn_MB(boost_matrix _genos_matrix, boost_matrix _phenos_matrix, int _subset_size, size_t _n_it_n, double _alpha_stat, list<unsigned> mem_a/*, P*/)
{
    list<unsigned> _markov_blanket_a;
    bool _markov_blanket_modified = true;
    int j = 0;
    forward(_markov_blanket_modified, _markov_blanket_a, _n_it_n, j/*, P*/, _genos_matrix, _phenos_matrix, _subset_size );
    return _markov_blanket_a;
}
//=================================================
// smmb_aco : forward
//=================================================
void smmb_aco::forward(bool _markov_blanket_modified, list<unsigned> _markov_blanket_a, size_t _n_it_n, int j/*, P*/, boost_matrix _genos_matrix, boost_matrix _phenos_matrix, int _subset_size )
{
    while (_markov_blanket_modified || (_markov_blanket_a.empty() && j<_n_it_n))
    {
        _markov_blanket_modified = false;
        /*
        TODO
        S = echantillone(P, _genos_matrix, k)
        s<-arg_max{score_association(s',T,_markov_blanket_a,mem_a)} //l'argument qui maximise
        */
        //if (p_valeur(s) << _alpha_stat) //TODO: Il faut une fonction pour calculer/renvoyer la p_valeur de la solution
        //{//cas de rejet de H_0
            //_markov_blanket_a = std::set_union (_markov_blanket_a, _markov_blanket_a + _markov_blanket_a.size(), S, S+S.size(), _markov_blanket_a.begin());// union de MB_a et S et retourne avec v.begin le debut du vecteur MB_a
            _markov_blanket_modified = true;
            backward(_markov_blanket_a, _phenos_matrix, _alpha_stat);
        //}
        j++;
    }
}
//=================================================
// smmb_aco : backward
//=================================================
void smmb_aco::backward(list<unsigned> _markov_blanket, boost_matrix _phenos_matrix, double _alpha_stat)
{
    for (size_t X = 0; X < _markov_blanket.size(); X++) {
        //for (size_t S = 0; S < count; S++) {
        //TODO: pour toute combinaison S non_vides inclus dans MB
            //independance_test_conditionnal(X,T,S_0); //TODO: omg c est chaud ca
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
void smmb_aco::run(boost_matrix _genos_matrix, boost_matrix _phenos_matrix, int _subset_size, size_t _n_it_n, size_t _n_ant, float _tau_0)
{
    // Initialization of Markov Blanket
    list<unsigned> _markov_blanket_s;
    float tau = _tau_0;
    list<unsigned> _markov_blanket_a;
    for (size_t i = 0; i < _n_it_n; i++)
    {
        //P <- calculer distribution, probabilité (tau, eta, _alpha_stat, beta)
        // For every ants
        for (size_t a = 0; a < _n_ant; a++)
        { // a parallelise
            // _genos_matrix<- echantillonner(P, D, KI) //TODO
            // Initialization of memory
            list<unsigned> mem_a;
            // Generate Markov Blanket and stock it in a temp variable
            _markov_blanket_a = learn_MB(_genos_matrix, _phenos_matrix, _subset_size, _n_it_n, _alpha_stat, mem_a/*, P*/);
        }
        // Initialization of final variable
        list<unsigned> mem;
        for (size_t a = 0; a < _n_ant; a++) {
            //ajouter(mem, mem_a); //TODO
            if (!_markov_blanket_a.empty()) {
                _markov_blanket_s.splice(_markov_blanket_s.end(), _markov_blanket_a); // move at the end of MB_s MB_a alternatively we can do MB_s.insert(MB_s.end(), MB_a.begin(), MB_a.end()); to copy MB_a content at MB_s end
            }
            //post traitement; //TODO
        }
    }
}
