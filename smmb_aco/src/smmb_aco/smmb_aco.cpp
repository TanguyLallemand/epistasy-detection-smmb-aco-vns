/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "smmb_aco.hpp"

//==============================================================================
// smmb_aco : constructor
//==============================================================================
smmb_aco::smmb_aco(data_parsing dataset, parameters_parsing params)
{
    // Storing datas in members objects
    this->_genos_matrix = dataset._geno_matrix;
    this->_pheno_vector = dataset._pheno_vector;
    this->_snp_id = dataset._snp_id_vector;
    this->_filename = dataset._geno_filename;
    // Miscellaneous
    this->_verbose = params.verbose;
    // Storing parameters in members variables
    this->_n_it_n = params.aco_n_iterations;
    this->_n_it = params.n_trials_to_learn_1_mb;
    this->_n_ant = params.aco_n_ants;
    this->_rho = params.aco_rho;
    this->_lambda = params.aco_lambda;
    this->_alpha_phero = params.aco_alpha;
    this->_beta_phero = params.aco_beta;
    this->_alpha_stat = params.alpha;
    this->_subset_size = params.aco_set_size;
    this->_sub_subset_size = params.subset_size_small;
    this->_output_prefix = params.output_prefix;
    this->_output_directory = params.output_directory;
    this->_pass_number = params.n_smmb_aco_runs;
    this->_tau_0 = params.aco_tau_init;

    // Initialization of the rng seed
    this->_rng.seed(time(NULL));

    // Initialization of vectors for pheromons
    this->_eta = boost_vector_float(_genos_matrix.size2(), params.aco_eta);
    this->_tau = boost_vector_float(_genos_matrix.size2(), _tau_0);
    this->_pheromone_distrib = boost_vector_float(_genos_matrix.size2(), 0.0);

    // Initialisation of the _markov_blanket_a to number of ant
    this->_markov_blanket_a.resize(_n_ant, false);

    // Initialisation of the _mem_ant to number of ant
    this->_mem_ant.resize(_n_ant, false);

    // Initialization of the distribution for the first iteration
    update_pheromon_distrib();
}

//==============================================================================
// smmb_aco : run
//==============================================================================
void smmb_aco::run()
{
    // Get actual time to measure running time
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    if (_verbose)
    {
        print_parameters();
    }
    std::cout << "### Backtrace of SMMB_ACO run ###" << '\n';
    // Loop, first level
    for (size_t i = 0; i < _n_it_n; i++)
    {
        std::cout << "Iteration #" << i << '\n';
        // Some optional verbose
        if (_verbose)
        {
            std::cout << "Tau vector" << '\n';
            std::cout << _tau << '\n';
        }
        // Reseting _markov_blanket_a for the new iteration
        _markov_blanket_a.clear();
        // Parallelization command, need to add correct argument in compilation line
        #pragma omp parallel for
        // Iterating through ants
        for (size_t a = 0; a < _n_ant; a++)
        {
            // Container for the ant subset
            boost_vector_int ant_subset(_subset_size);
            // Assigning a subset of _subset_size SNPs
            ant_subset = tools::sampling(_subset_size, _pheromone_distrib, _rng);
            // Generate MB from the ant subset
            learn_MB(ant_subset, _markov_blanket_a(a), _mem_ant(a));
        }
        save_iteration_result();
        // Update pheromon based on the results of tests performed on this iteration
        update_tau();
    }
    // prepare the second pass
    // listing all different SNPs in MB_s
    list<unsigned> new_set;
    // Iterate tought markov blanket
    for (auto d : _markov_blanket_s)
    {
        for (auto d2 : d.first)
        {
            // If a SNP is found in markov blanket and is not already in new dataset
            if (find(new_set.begin(), new_set.end(), d2) != new_set.end())
            {
                // Add it in data set
                new_set.push_back(d2);
            }
        }
    }
    // If user asked for two pass and if dataset is enought wide
    if ((_pass_number > 1) && (new_set.size() > 20))
    {
        next_pass(new_set);
    }
    // Final part of the algorithm, producing result file
    else
    {
        // Calculate scores of patterns in _markov_blanket_s into _stats_results
        score_for_final_results();
        // In order to print results of the run in terminal, please uncomment next line
        // show_results();
        // Save time
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        // Calculate time of execution
        this -> _duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        // Post treatment
        save_results();
    }
}

//==============================================================================
// smmb_aco : update_tau
//==============================================================================
void smmb_aco::update_tau()
{
    // Evaporate all the SNPs
    for (size_t i = 0; i < _tau.size(); i++)
    {
        _tau(i) = (1-_rho) * _tau(i);
    }
    // Iterate through memory map
    for (auto const& it : _mem)
    {
        for (auto const& score : it.second)
        {
            // Add pheromon based on the score for each SNP in _mem
            _tau(it.first) = _tau(it.first) + (_lambda * score);
        }
    }
    // Repercuting changes on pheromon distribution
    update_pheromon_distrib();
}

//==============================================================================
// smmb_aco : update_pheromon_distrib
//==============================================================================
void smmb_aco::update_pheromon_distrib()
{
    // Iterate tought pheromon distribution and update every values
    for (size_t i = 0; i < _pheromone_distrib.size(); i++)
    {
        _pheromone_distrib(i) = pow(_tau(i), _alpha_phero) + pow(_eta(i), _beta_phero);
    }
}

//==============================================================================
// smmb_aco : learn_MB
//==============================================================================
// Return non optimal Markov Blanket eventually empty
void smmb_aco::learn_MB(boost_vector_int & ant_subset, list<unsigned> & MB_a_ref, map<unsigned, list<float>> & mem_ant_ref)
{
    // Clearing memory of the predecent iteration
    mem_ant_ref.clear();
    // Initialization at true to enter the loop on first iteration
    bool markov_blanket_modified = true;
    // Counter of iterations
    unsigned j = 0;
    list<unsigned> save_MB;
    // Loop to generate the markov blanket
    while ((MB_a_ref.empty() && j<_n_it) && markov_blanket_modified)
    {
        // Saving MB to check differencies and launch next iteration or not
        save_MB = MB_a_ref;
        forward(MB_a_ref, ant_subset, mem_ant_ref);
        j++;
        // Sort markov blanket
        MB_a_ref.sort();
        // Update state of markov blanket
        if (save_MB != MB_a_ref)
        {
            markov_blanket_modified = true;
        }
        else
        {
            markov_blanket_modified = false;
        }
    }
}

//==============================================================================
// smmb_aco : forward
//==============================================================================
void smmb_aco::forward(list<unsigned> & MB_a_ref, boost_vector_int const& ant_subset, std::map<unsigned, list<float>> & mem_ant_ref)
{
    // Initialization of sub_subset container at right size with zeros
    boost_vector_int sub_subset(_sub_subset_size, 0);
    // Sub_sampling from ant_subset not taking SNPs already in the MB of this ant
    sub_sampling(sub_subset, ant_subset, MB_a_ref);
    // Container for all combinations of the current sub_subset
    list<list<unsigned>> combi_list;
    // Generating all the combinations from the drawn sub_subset
    get_all_combinations(sub_subset, combi_list);
    // Container to save the best pattern found by best_combination()
    list<unsigned> best_pattern;
    // Searching for the best combination based on score and returning the score and p_value of the solution
    boost_vector_float result = best_combination(best_pattern, combi_list, MB_a_ref, mem_ant_ref);
    // If independance hypothesis is rejected : entering here
    if (result(1) < _alpha_stat)
    {
        // Append the best pattern found at the end of the MB
        for (auto const& i : best_pattern)
        {
            MB_a_ref.push_back(i);
        }
        // Entering backward phase to remove worst SNPs of the MB
        backward(MB_a_ref);
    }
}

//==============================================================================
// smmb_aco : backward
//==============================================================================
void smmb_aco::backward(list<unsigned> & MB_a_ref)
{
    // Creating a copy of the MB to avoid modifying the iterator
    list<unsigned> iterate = MB_a_ref;
    // iterating SNPs of the MB
    for (auto const& current_SNP : iterate)
    {
        // Creating a copy of the current MB
        list<unsigned> MB_minus_current_SNP = MB_a_ref;
        // Removing the current SNP from the MB copy
        MB_minus_current_SNP.remove(current_SNP);
        // Converting current MB without current SNP to a vector
        boost_vector_int vector_MB(MB_minus_current_SNP.size());
        // Initialize a counter for index to copy MB into a vector
        int i=0;
        // Stocking MB in vector
        for (auto const& value : MB_minus_current_SNP)
        {
            vector_MB(i) = value;
            i++;
        }
        // Stocking the list of combination for the MB minus current SNP
        list<list<unsigned>> combi_list;
        // Generating all the combination from the drawn sub_subset
        get_all_combinations(vector_MB, combi_list);
        //reference to the column of the current SNP
        boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);
        //iterating the combination list
        for (auto const& combi : combi_list)
        {
            // Return an array with in first cell, g2 score and in second cell, associated p_value
            boost_vector_float result = statistics::make_contingencies_g_2_conditional_test_indep(mc, _pheno_vector, combi);
            // If the p-value is not significant entering here
            if (result(1) > _alpha_stat)
            {
                // Remove current SNP from the MB
                MB_a_ref.remove(current_SNP);
                break;
            }
        }
    }
}


//==============================================================================
// smmb_aco : sub_sampling
//==============================================================================
void smmb_aco::sub_sampling(boost_vector_int & sub_subset, boost_vector_int const& ant_subset, list<unsigned> MB_a_ref)
{
    //sub weight vector associated to the ant_subset
    boost_vector_float small_distrib(_subset_size);
    // Getting _tau values associated to the ant_subset SNPs
    for (size_t i = 0; i < ant_subset.size(); i++)
    {
        //if SNPs are already in the ant MB we won't pick them again
        if ((std::find(MB_a_ref.begin(), MB_a_ref.end(), ant_subset(i)) != MB_a_ref.end()))
        {
            small_distrib(i) = 0;
        }
        else
        {
            small_distrib(i) = _pheromone_distrib(ant_subset(i));
        }
    }
    // Picking SNPs in the ant subset
    boost_vector_int temporary = tools::sampling(_sub_subset_size, small_distrib, _rng);
    // Taking the SNPs on index returned by tools::sampling in ant_subset
    for (size_t j = 0; j < _sub_subset_size; j++)
    {
        sub_subset (j) = ant_subset(temporary(j));
    }
}

//==============================================================================
// smmb_aco : get_all_combinations
//==============================================================================
void smmb_aco::get_all_combinations(boost_vector_int & sub_subset, list<list<unsigned>> & combi_list)
{
    // Convert vector into list
    list<unsigned> subset(sub_subset.begin(), sub_subset.end());
    // Temporary list that we will append to combi list every time it is modified
    list<unsigned> temp;
    // Recursive function to generate all non-empty combinations
    generate_combinations(temp, combi_list, subset);
}

//==============================================================================
// smmb_aco : generate_combinations
//==============================================================================
// Recursively generate all non empty possible combinations from a set of snp
void smmb_aco::generate_combinations(list<unsigned> temp, list<list<unsigned>> & combi_list, list<unsigned> subset)
{
    // Copy the subset for next iteration
    list<unsigned> next_subset(subset);
    // Iterate the subset list
    for (auto const& h : subset)
    {
        // Add current SNP to the temp list
        temp.push_back(h);
        // Stocking the temp in the list of combinations
        combi_list.push_back(temp);
        // Remove the current x
        next_subset.remove(h);
        // Recursive call on the list without current x
        generate_combinations(temp, combi_list, next_subset);
        // Remove predecent SNP
        temp.pop_back();
    }
}

//==============================================================================
// smmb_aco : best_combination
//==============================================================================
boost_vector_float smmb_aco::best_combination(list<unsigned> & best_pattern, list<list<unsigned>> const& pattern_list, list<unsigned> & MB_a_ref, std::map<unsigned, list<float>> & mem_ant_ref)
{
    // Stock the current best_result
    boost_vector_float best_result(3, 1);
    // Iterate through the list of pattern
    for (auto const& current_pattern : pattern_list)
    {
        boost_vector_float result_pattern(3, 1);
        for (auto const& current_SNP : current_pattern)
        {
            // Setting up the list of conditionnals SNPs
            list<unsigned> conditionnal_set = current_pattern;
            MB_a_ref.sort();
            conditionnal_set.sort();
            conditionnal_set.merge(MB_a_ref);
            conditionnal_set.remove(current_SNP);
            // Making a matrix column ref to pass to the function
            boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, current_SNP);
            // Calculating score of the current SNP of the pattern
            boost_vector_float result_SNP = statistics::make_contingencies_g_2_conditional_test_indep(mc, _pheno_vector, conditionnal_set);
            // Stocking result in the ant_memory
            mem_ant_ref[current_SNP].push_back(result_SNP(0));
            // Score of the pattern
            if (result_SNP(1) < result_pattern(1))
            {
                result_pattern = result_SNP;
            }
        }
        // If the score of the current pattern is better than the current best it become the best pattern
        if (result_pattern(1) < best_result(1))
        {
            best_result = result_pattern;
            best_pattern = current_pattern;
        }
    }
    return best_result;
}

//==============================================================================
// smmb_aco : save_iteration_result
//==============================================================================
// save the result of each ant for this iteration in the global results
void smmb_aco::save_iteration_result()
{
    // Clearing the global memory
    _mem.clear();
    // Initialization of final variables of the run (merging all ants result)
    for (size_t a = 0; a < _n_ant; a++)
    {
        // Merging memory of all ants into global memory
        for (auto const& it : _mem_ant(a))
        {
            // Iterate map of the current ant
            for (auto const& it2 : it.second)
            {
                // Save map of the current ant in the global memory
                _mem[it.first].push_back(it2);
            }
        }

        // If ant's markov blanket is not empty save it into _markov_blanket_s
        if (!_markov_blanket_a(a).empty())
        {
            // Sort markov blanket of current ant
            _markov_blanket_a(a).sort();
            // if the ant MB is already known this add 1 to occurences number else it create the entry in the map and set it to 1
            _markov_blanket_s[_markov_blanket_a(a)] += 1;

        }
    }
}
//==============================================================================
// smmb_aco : next_pass
//==============================================================================
// prepare next pass by making unpickable non choosen SNPs and reinitialising other ones pheromon values
void smmb_aco::next_pass(list<unsigned> new_set)
{
    // Put n_smmb_aco_runs to 1 to stop above next time
    _pass_number--;
    for (auto v : new_set)
    {
        // Search for SNPs contained in new set
        if (find(new_set.begin(), new_set.end(), v) != new_set.end())
        {
            // If SNPs is found reinit his tau value to 100
            _tau[v] = _tau_0;
        }
        else
        {
            // If SNPs are not in dataset used for the second pass of SMMB_ACO their pheromons are set to 0 : they can't be picked anymore
            _tau[v] = 0;
            _eta[v] = 0;
        }
    }
    // Reset residual markov blanket
    _markov_blanket_s.clear();

    update_pheromon_distrib();
    // Run a new SMMB-ACO using first SMMB-ACO's output as input
    run();
}
//==============================================================================
// smmb_aco : score_for_final_results
//==============================================================================
// Calculating the score of the final selected patterns to print them in output file
void smmb_aco::score_for_final_results()
{
    _stats_results.resize(_markov_blanket_s.size());
    unsigned st=0;
    for (auto pattern : _markov_blanket_s)
    {
        _stats_results(st) = boost_vector_float(3, 1);
        for (auto snp : pattern.first)
        {
            boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, snp);
            list<unsigned> conditionnal_set = pattern.first;
            conditionnal_set.remove(snp);
            boost_vector_float temp_res =  statistics::make_contingencies_g_2_conditional_test_indep(mc, _pheno_vector, conditionnal_set);
            if (temp_res(1) < _stats_results(st)(1))
            {
                _stats_results(st) = temp_res;
            }
        }
        st++;
    }
}

//==============================================================================
// smmb_aco : show_results
//==============================================================================
void smmb_aco::show_results(vector<pair<unsigned, float>> & sorting, vector<string> & string_pat)
{
    cout << endl;
    cout << "# Result from SMMB-ACO \n";
    cout << "# Pattern || Occurences || G2-score || p-value || unreliable case\n";
    for (auto ss : sorting)
    {
        std::cout << string_pat[ss.first];
    }
    cout << endl;
}



//==============================================================================
// smmb_aco : save_results
//==============================================================================
void smmb_aco::save_results()
{
    // Parse filename
    // search for last "/"
    size_t firstindex = _filename.find_last_of("/");
    string filename_without_extension = _filename.substr(firstindex+1, 545754);
    // Search for "." and get index
    size_t lastindex = filename_without_extension.find_last_of(".");
    // Remove extension of given filename
    filename_without_extension = filename_without_extension.substr(0, lastindex);
    // Create the output file
    ofstream output_file(_output_directory + _output_prefix + filename_without_extension + "_smmb_aco.txt");
    // Fill output file with pattern, number of occurence and associated g2 score and p-value
    output_file << "# Result from SMMB-ACO \n";
    output_file << "# Pattern || Occurences || G2-score || p-value || unreliable case\n";

    vector<pair<unsigned, float>> sorting;
    vector<string> string_pat;

    format_results(sorting, string_pat);

    for (auto ss : sorting)
    {
        output_file << string_pat[ss.first];
    }
    // Print in file time of exectution
    output_file << "# Time of execution: " << _duration << " milliseconds" << endl;
    // Verbose displayed in terminal
    if (_verbose)
    {
        show_results(sorting, string_pat);
    }
    std::cout << "# Time of execution: " << _duration << " milliseconds" << endl;
    std::cout << "### SMMB_ACO has finished please see results in: " << '\n' << _output_directory + _output_prefix + filename_without_extension + "_smmb_aco.txt" << '\n';
}

//==============================================================================
// smmb_aco : compareFunc
//==============================================================================
bool smmb_aco::compareFunc(pair<unsigned, float> const& a, pair<unsigned, float> const& b)
{
    if (a.second < b.second)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//==============================================================================
// smmb_aco : print_parameters
//==============================================================================
void smmb_aco::print_parameters()
{
    std::cout << "### Parameters used for this run: " << endl;
    std::cout << "Number of iterations in ACO: " << _n_it_n  << endl;
    std::cout << "Maximal number of iterations allowed to learn one Markov blanket: " << _n_it  << endl;
    std::cout << "Number of ants: " << _n_ant  << endl;
    std::cout << "Evaporation rate: " << _rho  << endl;
    std::cout << "Lambda parameter from ACO-PDF update function: " << _lambda  << endl;
    std::cout << "Alpha, used to ajust the relative importance between pheromone rate and a priori knowledge: " << _alpha_phero  << endl;
    std::cout << "Beta, used to ajust the relative importance between pheromone rate and a priori knowledge: " << _beta_phero  << endl;
    std::cout << "Alpha type I error rate: " << _alpha_stat  << endl;
    std::cout << "Number of snps sampled in each ant: " << _subset_size  << endl;
    std::cout << "Size of the smallest subset: " << _sub_subset_size  << endl;
    std::cout << "Prefix of output files: " << _output_prefix  << endl;
    std::cout << "Path of output directory: " << _output_directory  << endl;
    std::cout << "Number of consecutive runs of SMMB-ACO: " << _pass_number  << endl;
    std::cout << "Value to initiate evaporation rates: " << _tau_0  << endl;
    std::cout << "### End of parameters" << endl;
    std::cout << endl;
}


//==============================================================================
// smmb_aco : format_results
//==============================================================================
void smmb_aco::format_results(vector<pair<unsigned, float>> & sorting, vector<string> & string_pat)
{
    unsigned tu = 0;
    unsigned to = 0;
    string lign;
    for (auto const& pattern : _markov_blanket_s)
    {
        // Only output results considered as reliable following statistical treshold fixed by user
        if (_stats_results(tu)[1] < _alpha_stat)
        {
            // Construct using concatenation output
            lign = "";
            lign = lign + '{';
            for (auto const& snp : pattern.first)
            {
                lign += _snp_id(snp);
                if (snp!=pattern.first.back()) {
                    lign = lign + ',';
                }
            }
            lign = lign + "} || ";
            lign = lign + to_string(pattern.second);
            lign = lign + " || ";
            // Add to output g2 score
            lign = lign + to_string(_stats_results(tu)[0]);
            lign = lign + " || ";
            stringstream ss;
            // Add to output p value using scientific writing and two significant digits
            ss << scientific << setprecision(2) << _stats_results(tu)[1];
            lign = lign + ss.str();
            lign = lign + " || ";
            // Add to output number of unreliable g2 tests
            lign = lign + to_string((int)_stats_results(tu)[2]);
            lign = lign + "\n";
            string_pat.push_back(lign);
            // Sort result following p-value in descending order
            pair<unsigned, float> pair_sort(to, _stats_results(tu)[1]);
            sorting.push_back(pair_sort);
            to++;
        }
        tu++;
    }
    sort(sorting.begin(), sorting.end(), compareFunc);
}
