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
void vns::neighborhood_change(list<unsigned> x, list<unsigned> second_x, int k)
{

}

//==============================================================================
//vns : shake
//==============================================================================
list<unsigned> vns::shake(list<unsigned> x)
{
    int index_pattern = rand() % (_neighborhood[x][0].size());

    list<unsigned> x2 = _neighborhood[x][0][index_pattern];

    return x2;
}

//==============================================================================
//vns : variable_neighborhood_descent
//==============================================================================
float vns::variable_neighborhood_descent(list<unsigned> second_x, list<unsigned> & third_x) // This is the VND phase
{
    float score, best_score = 0;

    for (auto neighbor_iterator : _neighborhood[second_x][0])
    {
        //TODO calcul du score de neighbor_iterator

        if (score > best_score)
        {
            best_score = score;
            third_x = neighbor_iterator;
        }
    }
    return best_score;
}

//==============================================================================
//vns : run
//==============================================================================
void vns::run()
{
    //initialisation of patterns and their neighbors
    generate_patterns();

    for (size_t i = 0; i < _n_it_max; i++)
    {
        //selecting starting pattern
        int index_pattern = rand() % (pattern_list.size());

        //initialisation of starting pattern
        list<unsigned> x = pattern_list[index_pattern];
        //TODO ici on doit tester le pattern de départ
        float x_score = 0;

        list<unsigned> second_x;
        list<unsigned> third_x;

        float third_x_score = 0;

        int iterator = 0;

        int k = 0;
        while (k < _k_max)
        {
            //take a random neighbor of x
            second_x = shake(x);

            //searching for the best neighbor of second_x
            third_x_score = variable_neighborhood_descent(second_x, third_x);

            if (third_x_score > x_score)
            {
                x = third_x;
                x_score = third_x_score;
                k = 0;
            }
            else
            {
                k++;
            }
        }

        //saving the local optimum
        auto current_opti = _optimum_set.find(x);
        if (current_opti->second.size() == 0)
        {
            current_opti->second = {1,x_score};
        }
        else
        {
            current_opti->second[0] += 1;
            current_opti->second[1] = x_score;
        }

    }
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
    //if we are on the _k_max recursive call we don't go deeper
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

            //add the pattern to the vector
            pattern_list.push_back(temp);

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
