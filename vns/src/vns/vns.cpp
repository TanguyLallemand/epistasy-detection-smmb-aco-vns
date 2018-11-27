#include "vns.hpp"

vns::vns(data_parsing dataset, parameters_parsing _params)
{
    //params unpacking
    this->_n_it_max = _params._n_it_max;
    this->_k_max = _params._k_max;

    //unpacking datas
    this->_genos_matrix = dataset._geno_matrix;
    this->_pheno_vector = dataset._pheno_vector;
    this->_snp_id = dataset._snp_id_vector;
    this->_filename = dataset._geno_filename.substr(14, dataset._geno_filename.length());
}

//==============================================================================
//vns : neighborhood_change
//==============================================================================
void vns::neighborhood_change(int x, int second_x, int k)
{

}

//==============================================================================
//vns : shake
//==============================================================================
void vns::shake(int x, int k)
{

}

//==============================================================================
//vns : variable_neighborhood_descent
//==============================================================================
void vns::variable_neighborhood_descent(int x, int k_max) // This is the VND phase
{

}

//==============================================================================
//vns : run
//==============================================================================
void vns::run(int x, int l_max, int k_max, int n_it_max)
{
    //initialisation of patterns and their neighbors
    generate_patterns();

    //selecting starting pattern
    int testi = rand() % (_neighborhood.size()-1);



    int iterator = 0;
    int second_x;
    int third_x;
    while (_n_it_max > iterator)
    {
        int k = 1;
        while (k != _k_max)
        {
            shake(x, k);
            variable_neighborhood_descent(second_x, l_max);
            neighborhood_change(x, third_x, k);
        }
        iterator++;
    }
    //return x;
}

//==============================================================================
//vns : generate_patterns
//==============================================================================
void vns::generate_patterns()
{
    list<unsigned> temp;
    list<unsigned> snp_list(_genos_matrix.size2());
    iota(snp_list.begin(), snp_list.end(), 0);
    generate_patterns(temp, snp_list);
    set_neighbors();
}

//==============================================================================
//vns : generate_patterns
//==============================================================================
void vns::generate_patterns(list<unsigned> temp, list<unsigned> snp_list)
{
    //TODO changer le 3 par une variable globale qu'on passe en arg pour la taille max du pattern recherché
    //if we are on the size_pattern recursive call we don't go deeper
    if (temp.size() < _k_max)
    {
        for (auto snp : snp_list)
        {
            //add the snp to the pattern
            temp.push_back(snp);

            //create the entry in the map for current pattern
            _neighborhood[temp];

            //copy the list of snp
            list<unsigned> next_snp_list = snp_list;

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

//==============================================================================
//vns : set_neighbors
//==============================================================================
void vns::set_neighbors()
{
    //iterating all patterns
    for (auto current : _neighborhood)
    {
        //iterate to generate list of pattern on a certain distance
        for (size_t i = 0; i < _k_max; i++)
        {
            //each box of the vector contains neighbor with a different distance to the associated pattern
            current.second.resize(_k_max);
            //iterating all patterns to find neighbors for current at i+1 distance
            for (auto candidat_neighbor : _neighborhood)
            {
                //finding the number of common values needed to be neighbors at current distance
                unsigned neighbors_score = max(current.first.size(), candidat_neighbor.first.size())-i-1;

                //counting common values between vectors
                unsigned common_values = 0;
                for (auto s : current.first)
                {
                    auto test = find(candidat_neighbor.first.begin(), candidat_neighbor.first.end(), s);
                    if (*test != s)
                    {
                        common_values+=1;
                    }
                }

                //if the 2 sets have only 1 difference they are neighbors
                if (neighbors_score == common_values)
                {
                    //append the new neighbor to the list of neighbors with i+1 distance
                    current.second[i].push_back(candidat_neighbor.first);
                }
            }
        }
    }
}
