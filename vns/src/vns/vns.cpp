#include "vns.hpp"

vns::vns(parameters_parsing _params)
{
    this->_n_it_max = _params.n_it_max;
    this->_k_max = _params.k_max;
}

void vns::neighborhood_change(int x, int second_x, int k)
{

}

void vns::shake(int x, int k)
{

}

void vns::variable_neighborhood_descent(int x, int k_max) // This is the VND phase
{

}

void vns::run(int x, int l_max, int k_max, int n_it_max)
{
    int iterator = 0;
    int second_x;
    int third_x;
    while (n_it_max > iterator)
    {
        int k = 1;
        while (k != k_max)
        {
            shake(x, k);
            variable_neighborhood_descent(second_x, l_max);
            neighborhood_change(x, third_x, k);
        }
        iterator++;
    }
    //return x;
}

void vns::generate_patterns()
{
    list<unsigned> temp;
    generate_patterns(temp, snp_list);
    set_neighbors();
}

void vns::generate_patterns(list<unsigned> temp, list<unsigned> snp_list)
{
    //TODO changer le 3 par une variable globale qu'on passe en arg pour la taille max du pattern recherché
    //if we are on the size_pattern recursive call we don't go deeper
    if (temp.size() < 3)
    {
        for (auto snp : snp_list)
        {
            //add the snp to the pattern
            temp.push_back(snp);

            //create the entry in the map for current pattern
            _neighborhood[temp];

            //copy the list of snp
            next_snp_list = snp_list;

            //remove snp added to temp pattern
            next_snp_list.remove(snp);

            //recursive call to get next pattern
            generate_patterns(temp, next_snp_list);

            //remove current snp from the list to let the next one come
            temp.pop_back();
        }
    }
    else
    {
        for (auto snp : snp_list)
        {
            //add the snp to the pattern
            temp.push_back(snp);

            //create the entry in the map for current pattern
            _neighborhood[temp];

            //remove current snp from the list to let the next one come
            temp.pop_back();
        }
    }
}

void vns::set_neighbors();
{
    //iterating all patterns
    for (auto current : _neighborhood)
    {
        //iterating all patterns to find neighbors for current
        for (auto candidat_neighbor : _neighborhood)
        {
            //finding the number of common values needed to be neighbors
            unsigned neighbors_score = max(current.first.size(), candidat_neighbor.first.size())-1;

            //counting common values between vectors
            unsigned common_values = 0;
            for (auto s : current.first)
            {
                if (find(candidat_neighbor.first.begin(), candidat_neighbor.first.end(), s))
                {
                    common_values+=1;
                }
            }

            //if the 2 sets have only 1 difference they are neighbors
            if (neighbors_score == common_values)
            {
                //append the new neighbor to the list
                current.second.push_back(candidat_neighbor.first);
            }
        }
    }
}
