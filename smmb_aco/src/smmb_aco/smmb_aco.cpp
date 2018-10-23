/*
D: Matrice de genotype
T: Vecteur de phénotype
n_it: nombre d'itération ACO
n_ants: nombre de fourmis
K: Taille du sous ensemble de variable echantilloné à partir de D, pour chaque fourmis
n_it_n: nombre d'itération maximale qui force l'exploration de l'espace de recherche dans Markov Blanket
alpha': seuil pour erreur de type I global
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
        MB_a <- learn_MB(D_a,T, k, N_it_n, alpha', mem_a, P)
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
        si(p_valeur(s)<alpha')//cas de rejet de H_0
        alors
            MB_a<-MB_a union s
            MB_modified <- true
            backward_phase(MB_a, T, aplpha')
        j++
    return MB_a
}

procedure backward(MB, T, alpha)
{
    pour tout x element de MB
        pour toute combinaison S non_vides inclus dans MB
            realiser le terst d'independance entre X et T conditionnellement a S_0
            Si p_valeur>alpha'//H_0: independance
            alors
            MB <- MB\{x} //a gerer en liste
            break
}
*/

smmb_aco()
{
    
}
