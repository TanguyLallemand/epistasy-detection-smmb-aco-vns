// A determiner: le type de x et de k
// Un peu de biblio

/*
C'est un vnd un peu ameliore via une perturbation

Selectionner une architecture de voisinnage N_k 1<=k<=k_max
creer une solution intiale
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
    retourner s, la meilleur solution rencontrer dans M
*/


// Definition de fonctions
/*
function generate_neighborhood()
{   //TODO detreminer le type de l espace des solution N_k
    1<=k<=k_max
    return N_k
}
function create_initial_solution
{   //TODO cette fonction consiste a tirer en random une solution au hasard

}
function neighborhood_change(x,x',k)
{
    if f(x') < f(x)
    {
        x = x' //Make a move
        k = 1 // Inital neighborhood
    }
    else
    {
        k = k+1 //Next neighborhood
    }
    return x,k
}
function Shake(x,k)
{
    w = (1+Rand(0,1) * N_k(s) )
    x' = x^w
    return x'
}

function VND(x, k_max)
{
    k = 1
    while (k != kmax)
    {
        x' = argmin f(y) //Find the best neighbor in N_k(x)
        x,k = neighborhood_change(x,x',k) //Change neighborhood
    }
    return x;
}

function general_VNS(x, l_max, k_max, n_it_max) //c'est le run
{
    while (n_it_max>it)
    {
        k = 1
        while (k != k_max)
        {
            x' = shake(x,k)
            x'' = VND(x',l_max)
            x,k = neighborhood_change(x,x'',k)
        }
        it++
    }
    return x;
}
*/
