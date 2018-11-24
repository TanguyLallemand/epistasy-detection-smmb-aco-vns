#include "smmb_aco.hpp"

//==============================================================================
// smmb_aco : constructeur
//==============================================================================
smmb_aco::smmb_aco(boost_matrix genos_matrix, boost_vector_int pheno_vector, parameters_parsing _params)
{
    //stocking datas in members objects
    this->_genos_matrix = genos_matrix;
    this->_pheno_vector = pheno_vector;

    //stocking parameters in members variables
    this->_n_it_n = _params.aco_n_iterations;
    this->_n_ant = _params.aco_n_ants;
    this->_rho = _params.aco_rho;
    this->_lambda = _params.aco_lambda;
    this->_alpha_phero = _params.aco_alpha;
    this->_beta_phero = _params.aco_beta;
    this->_alpha_stat = _params.alpha;
    this->_subset_size = _params.aco_set_size;
    this->_sub_subset_size = _params.subset_size_small;

    //Initialization of the rng seed
    this->_rng.seed(2);

    //Initialization of vectors for pheromons
    this->_eta = boost_vector_float(_genos_matrix.size2(), (float)_params.aco_eta);
    this->_tau = boost_vector_float(_genos_matrix.size2(), (float)_params.aco_tau_init);
    this->_pheromone_distrib = boost_vector_float(_genos_matrix.size2(), 0.0);

    //initialisation of the _markov_blanket_a to number of ant
    this->_markov_blanket_a.resize(_n_ant, false);

    //initialisation of the _mem_ant to number of ant
    this->_mem_ant.resize(_n_ant, false);

    //Initialization of the distribution for the first iteration
    update_pheromon_distrib();
}

//==============================================================================
// smmb_aco : update_tau
//==============================================================================
void smmb_aco::update_tau()
{
    //evaporate all the SNPs
    for (size_t i = 0; i < _tau.size(); i++) {
        _tau(i) = (1-_rho) * _tau(i);
    }

    //iterate through memory map
    for (auto const& it : _mem)
    {
        for (auto const& score : it.second) {
            //add pheromon based on the score for each SNP in _mem
            _tau(it.first) = _tau(it.first) + (_lambda * score);
        }
    }

    //repercuting changes on the distribution
    update_pheromon_distrib();
}

//==============================================================================
// smmb_aco : update_pheromon_distrib
//==============================================================================
void smmb_aco::update_pheromon_distrib()
{
    for (size_t i = 0; i < _pheromone_distrib.size(); i++) {
        _pheromone_distrib(i) = pow(_tau(i), _alpha_phero) + pow(_eta(i), _beta_phero);
    }
}
//==============================================================================
// smmb_aco : learn_MB
//==============================================================================
//Return non optimal Markov Blanket eventually empty
void smmb_aco::learn_MB(boost_vector_int & ant_subset, list<unsigned> & MB_a_ref, map<unsigned, list<float>> & mem_ant_ref)
{
    // std::cout << "learn_MB" << '\n';
    //clearing manually the ant memory map to avoid corruption of the object
    for (size_t i = 0; i < _tau.size(); i++) {
        mem_ant_ref[i].clear();
    }

    //to enter the loop on first iteration
    bool markov_blanket_modified = true;

    //counter of iterations
    unsigned j = 0;

    //loop to generate the markov blanket //TODO check the condition
    while (markov_blanket_modified || (!MB_a_ref.empty() && j<_n_it_n))
    {
        forward(markov_blanket_modified, MB_a_ref, ant_subset, mem_ant_ref);
        j++;
    }
    // std::cout << "fin learn_MB" << '\n';
    //backward(markov_blanket_modified, MB_a_ref);
    //TODO voir si on en met un la finalement
}

//==============================================================================
// smmb_aco : forward
//==============================================================================
void smmb_aco::forward(bool & markov_blanket_modified, list<unsigned> & MB_a_ref, boost_vector_int const& ant_subset, std::map<unsigned, list<float>> & mem_ant_ref)
{
    // std::cout << "forward" << '\n';
    //to break the loop if nothing modified
    markov_blanket_modified = false;

    // sub_subset container
    boost_vector_int sub_subset(_sub_subset_size, 0);

    // sub_sampling from ant_subset not taking SNPs already in the MB of this ant
    sub_sampling(sub_subset, ant_subset, MB_a_ref);

    //container for all combinations of the current sub_subset
    list<list<unsigned>> combi_list;

    // generating all the combination from the drawn sub_subset
    get_all_combinations(sub_subset, combi_list);

    //container to save the best pattern found by best_combination()
    list<unsigned> best_pattern;

    // searching for the best combination based on score and returning the score and p_value of the solution
    boost_vector_float result(2, 0.0);
    result = best_combination(best_pattern, combi_list, MB_a_ref, mem_ant_ref);
    //if independance hypothesis is rejected : entering here
    if (result(1) < _alpha_stat)
    {
        //append the best pattern found at the end of the MB
        for (auto const& i : best_pattern) {
            MB_a_ref.push_back(i);
        }
        //markov blanket has been modified we need an other loop run
        markov_blanket_modified = true;

        // entering backward phase to remove worst SNPs of the MB
        backward(MB_a_ref);
    }
    // std::cout << "fin forward" << '\n';
}

//==============================================================================
// smmb_aco : backward
//==============================================================================
void smmb_aco::backward(list<unsigned> & MB_a_ref)
{
    // std::cout << "backward" << '\n';
    //creting a copy of the MB to avoid modifying the iterator
    list<unsigned> iterate = MB_a_ref;
    // iterating SNPs of the MB
    for (auto const& current_SNP : iterate)
    {
        //creating a copy of the current MB
        list<unsigned> MB_minus_current_SNP = MB_a_ref;

        //removing the current SNP from the current MB
        MB_minus_current_SNP.remove(current_SNP);

        // converting current MB without current SNP to a vector
        boost_vector_int vector_MB(MB_minus_current_SNP.size());

        //counter for index to recopy MB into a vector
        int i=0;
        //stocking MB in vector
        for (auto const& value : MB_minus_current_SNP)
        {
            vector_MB(i) = value;
            i++;
        }
        //stocking the list of combination for the MB minus current SNP
        list<list<unsigned>> combi_list;

        //generating all the combination from the drawn sub_subset
        get_all_combinations(vector_MB, combi_list);

        //reference to the column of the current SNP
        boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);

        //iterating the combination list
        for (auto const& combi : combi_list)
        {
            //stocking result of the test
            boost_vector_float result(2, 0.0);

            // Return an array with in first cell, chi 2 score and in second cell, associated p_value
            result = statistics::make_contingencies_chi_2_conditional_test_indep(mc, _pheno_vector, combi);

            //if the p-value is not significant entering here
            if (result(1) > _alpha_stat)
            {
                //remove current SNP from the MB
                MB_a_ref.remove(current_SNP);
                break;
            }
        }
    }
    // std::cout << "fin backward" << '\n';
}

//==============================================================================
// smmb_aco : run
//==============================================================================
void smmb_aco::run()
{
    std::cout << "run" << '\n';

    for (size_t i = 0; i < _n_it_n; i++)
    {
        std::cout << "iteration # " << i << '\n';
        std::cout << "tau vector" << '\n';
        std::cout << _tau << '\n';

        //reinitialisation of the MB of each ant (avoiding core dump by doing it one by one)
        for (size_t k = 0; k < _markov_blanket_a.size(); k++) {
            _markov_blanket_a(i).clear();
        }

        //TODO For every ants a parallelise : #pragma omp parallel for
        //iterating through ants
        for (size_t a = 0; a < _n_ant; a++)
        {
            // std::cout << "ant counter" << '\n';
            // std::cout << a << '\n';
            //container for the ant subset
            boost_vector_int ant_subset(_subset_size);

            //assigning a subset of _subset_size SNPs
            ant_subset = tools::sampling(_subset_size, _pheromone_distrib, _rng);

            // std::cout << "ant_subset" << '\n';
            // std::cout << ant_subset << '\n';
            //generate MB from the ant subset
            learn_MB(ant_subset, _markov_blanket_a(a), _mem_ant(a));
        }

        // clearing the global memory (one by one to avoid core dump)
        for (size_t i = 0; i < _tau.size(); i++) {
            _mem[i].clear();
        }

        // Initialization of final variables
        for (size_t a = 0; a < _n_ant; a++)
        {
            //merging memory of all ants into global memory
            for (auto const& it : _mem_ant(a))
            {
                //iterate map of the current ant
                for (auto const& it2 : it.second)
                {
                    _mem[it.first].push_back(it2);
                }
            }

            // If ant's markov blanket is not empty
            if (!_markov_blanket_a(a).empty()) //FIXME
            {
                //save the ant MB ant at the end of MB list
                _markov_blanket_s.push_back(_markov_blanket_a(a));
            }
        }
        //add pheromon based on the score
        update_tau();
    }
    std::cout << "fin run" << '\n';
    //post treatment
}

//==============================================================================
// smmb_aco : sub_sampling
//==============================================================================
void smmb_aco::sub_sampling(boost_vector_int & sub_subset, boost_vector_int const& ant_subset, list<unsigned> MB_a_ref)
{
    // std::cout << "sub_sampling" << '\n';
    //sub weight vector associated to the ant_subset
    boost_vector_float small_distrib(_subset_size);

    //getting _tau values associated to the ant_subset SNPs
    for (size_t i = 0; i < ant_subset.size(); i++)
    {
        if ((std::find(MB_a_ref.begin(), MB_a_ref.end(), ant_subset(i)) != MB_a_ref.end()))
        {
            small_distrib(i) = 0;
        }
        else
        {
            small_distrib(i) = _pheromone_distrib(ant_subset(i));
        }
    }

    //stocking the position of picked SNP in ant_subset
    boost_vector_int temporary(_sub_subset_size, 0);

    //pisking SNPs in the ant subset
    temporary = tools::sampling(_sub_subset_size, small_distrib, _rng);

    //taking the SNPs on index returned by tools::sampling in ant_subset
    for (size_t j = 0; j < _sub_subset_size; j++)
    {
        sub_subset (j) = ant_subset(temporary(j));
    }
    // std::cout << "fin sub_sampling" << '\n';
}

//==============================================================================
// smmb_aco : get_all_combinations
//==============================================================================
void smmb_aco::get_all_combinations(boost_vector_int & sub_subset, list<list<unsigned>> & combi_list)
{
    // std::cout << "get_all_combinations" << '\n';
    //convert vector into list
    list<unsigned> subset(sub_subset.begin(), sub_subset.end());

    //temporary list that we will append to combi list every time it is modified
    list<unsigned> temp;

    //recursive function to generate all non-empty combinations
    generate_combinations(temp, combi_list, subset);
    // std::cout << "fin get_all_combinations" << '\n';
}

//==============================================================================
// smmb_aco : generate_combinations
//==============================================================================
void smmb_aco::generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset)
{
    // std::cout << "generate_combinations" << '\n';
    //copy the subset for next iteration
    list<unsigned> next_subset(subset);

    //iterate the subset list
    for (auto const& h : subset)
    {
        //add current snp to the temp list
        temp.push_back(h);
        //stocking the temp in the list of combinations
        combi_list.push_back(temp);
        //remove the current x
        next_subset.remove(h);
        //recursive call on the list without current x
        generate_combinations(temp, combi_list, next_subset);
        //remove predecent snp
        temp.pop_back();
    }
    // std::cout << "fin generate_combinations" << '\n';
}

//==============================================================================
// smmb_aco : best_combination
//==============================================================================
boost_vector_float smmb_aco::best_combination(list<unsigned> & best_pattern, list<list<unsigned>> const& pattern_list, list<unsigned> & MB_a_ref, std::map<unsigned, list<float>> & mem_ant_ref)
{
    // std::cout << "best_combination" << '\n';
    //stock the current best_result
    boost_vector_float best_result(2, 0);
    //iterate through the list of pattern
    for (auto const& current_pattern : pattern_list)
    {
        boost_vector_float result_pattern(2, 0);
        for (auto const& current_SNP : current_pattern)
        {
            //setting up the list of conditionnals SNPs
            list<unsigned> conditionnal_set = current_pattern;
            MB_a_ref.sort();
            conditionnal_set.sort();
            conditionnal_set.merge(MB_a_ref);
            conditionnal_set.remove(current_SNP);
            //making a matrix column ref to pass to the function

            boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);

            //calculating score of the current SNP of the pattern and add it to the pattern score
            boost_vector_float result_SNP(2, 0.0);
            result_SNP = statistics::make_contingencies_chi_2_conditional_test_indep(mc, _pheno_vector, conditionnal_set);
            //stocking result in the ant_memory
            mem_ant_ref[current_SNP].push_back(result_SNP(0));

            //score of the pattern
            if (result_SNP(0) > result_pattern(0))
            {
                result_pattern = result_SNP;
            }
        }
        //if the score of the current pattern is better than the current one it become the current best pattern
        if (result_pattern(0) > best_result(0))
        {
            best_result = result_pattern;
            best_pattern = current_pattern;
        }
    }
    // std::cout << "fin best_combination" << '\n';
    return best_result;
}
