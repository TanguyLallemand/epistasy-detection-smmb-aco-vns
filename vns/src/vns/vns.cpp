/*
C'est un vnd un peu ameliore via une perturbation

Selectionner une architecture de voisinnage N_k 1<=k<=k_max
creer une solution intiale scanf("initialisation int main(int argc, char const *argv[]) {
    tant que non condition fin
        k<-1
        tant que k<k_max
            choisis aleatoiremùent un voisin s' dans N_k(s)
            s''<-localsearch(s',N_k)
            si (score(s'')<score(s))
                s<-s'', k=1
            sinon
                k++
        memoriser l'optimum local obtenue dans M
    retourner s, la meilleur solution rencontrer dans M*/
