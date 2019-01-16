/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
#include "vns.hpp"

//==============================================================================
// vns : vns
//==============================================================================
//Initialization of all variable for the run

vns::vns(data_parsing dataset, parameters_parsing _params)
{
// Parameters unpacking
    // Output informations
    this->_output_directory = _params.output_directory;
    this->_output_prefix = _params.output_prefix;
    this->_verbose = _params.verbose;
    // Algorithm parameters
    this->_iteration_num = _params._iteration_num;
    this->_alpha = _params.alpha;
    this->_pat_size_max = _params._pat_size_max;
    this->_pat_size_min = _params._pat_size_min;
    this->_max_it_vns = _params._max_it_vns;
    this->_max_it_local_search = _params._max_it_local_search;
    // Definition of neighbor exploration range
    this->_k_max = _params._pat_size_max-1;
    this->_l_max = _params._pat_size_max-1;

// Unpacking datas
    this->_genos_matrix = dataset._geno_matrix;
    this->_phenos_vector = dataset._pheno_vector;
    this->_snp_id = dataset._snp_id_vector;
    this->_filename = dataset._geno_filename;

// Initialization of the rng seed for all random picks
    this->_rng.seed(time(NULL));
}

//==============================================================================
// vns : run
//==============================================================================

void vns::run()
{
    // Initialization of time
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    print_parameters();
    std::cout << "### Backtrace of VNS run ###" << '\n';
    // Parallelization
    #pragma omp parallel for
    for (size_t i = 0; i < _iteration_num; i++)
    {
        std::cout << "Iteration # : " << i << '\n';
        // Selecting starting pattern
        vector<unsigned> x = generate_starting_pattern();
        // Compute the score of x for first iteration
        vector<float> x_score = test_pattern(x);
        // Declaration of temporary storing
        vector<unsigned> second_x, third_x;
        vector<float> third_x_score;
        // Stop when x did not changed for k_max itérations
        unsigned exploration = 0;
        int k = 1;
        // Exploring neighborhood to a range of _k_max
        while (k < _k_max)
        {
            // Exploration of the current neighborhood
            while (exploration < _max_it_vns)
            {
                // Take a random neighbor of x
                second_x = shake(x, k);
                // Searching for a best neighbor of second_x
                third_x_score = local_search(second_x, third_x);
                if (third_x_score[1] < x_score[1])
                {
                    x = third_x;
                    x_score = third_x_score;
                    k = 1;
                    exploration = 0;
                }
                else
                {
                    exploration++;
                }
            }
            save_local_optimum(x, x_score);
            k++;
            exploration = 0;
        }
        // Saving the local optimum
    }
    // End time
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    this-> _duration = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    // Write results in a file
    write_result_file();
}


//==============================================================================
// vns : generate_starting_pattern
//==============================================================================
//Take a pattern with a size between _pat_size_min and _pat_size_max and fill it
// With random SNPs

vector<unsigned> vns::generate_starting_pattern()
{
    // Initialization of list to store the pattern
    vector<unsigned> pattern;
    // Random pick of the pattern size between _pat_size_min and _pat_size_max
    std::uniform_int_distribution<int> distribution(_pat_size_min,_pat_size_max);
    unsigned size_pattern = distribution(_rng);
    // Fill the pattern with size_pattern different SNPs
    for (size_t i = 0; i < size_pattern; i++)
    {
        // Pick a random SNP
        std::uniform_int_distribution<int> distribution(0,_genos_matrix.size2()-1);
        unsigned new_SNP = distribution(_rng);
        // Make sure to pick unique SNPs
        while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
        {
            new_SNP = distribution(_rng);
        }
        // Add the selected SNP to the pattern
        pattern.push_back(new_SNP);
    }
    return pattern;
}

//==============================================================================
// vns : local_search
//==============================================================================
// Exploration of the neighborhood of the provided pattern
vector<float> vns::local_search(vector<unsigned> second_x, vector<unsigned> & third_x)
{
    // Initialization of temporary variables
    vector<float> score;
    vector<float> best_score = test_pattern(second_x);
    third_x = second_x;
    unsigned l=1;
    vector<unsigned> candidat_neighbor;
    unsigned exploration = 0;
    while (l < _l_max)
    {
        while (exploration < _max_it_local_search)
        {
            candidat_neighbor = shake(third_x, l);
            score = test_pattern(candidat_neighbor);
            if (score[1] < best_score[1])
            {
                // If this pattern is the best tested he becomes the new best and reset counters
                best_score = score;
                third_x = candidat_neighbor;
                l=1;
                exploration = 0;
            }
            else
            {
                // If this pattern isn't better than the current we keep the current and increment the counter
                exploration++;
            }
        }
        // Get to the next neighborhood
        l++;
        // Reset the neightborhood exploration counter for next scale neighborhood
        exploration = 0;
    }
    return best_score;
}

//==============================================================================
// vns : test_pattern
//==============================================================================
//test the pattern and return a vector containing score,p-value,unreliable case

vector<float> vns::test_pattern(vector<unsigned> const& pattern)
{
    // Prepare the list of matrix column for the test
    vector<boost::numeric::ublas::matrix_column<boost_matrix>> pattern_datas;
    for (auto snp : pattern)
    {
        boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, snp);
        pattern_datas.push_back(mc);
    }
    // Test the datas provided
    vector<float> result = statistics::compute_p_value(pattern_datas, _phenos_vector);
    return result;
}

//==============================================================================
// vns : shake
//==============================================================================
// Return a neighbor of pattern, this neighbor have a distance of k with
// the pattern

vector<unsigned> vns::shake(vector<unsigned> pattern, unsigned k)
{
    unsigned mutation_type;
    uniform_int_distribution<int> distribution_pattern(0,pattern.size()-1);
    uniform_int_distribution<int> distribution_snp(0,_genos_matrix.size2()-1);
    // Make sure we won't do a forbiden change in the pattern
    if (_pat_size_min == _pat_size_max)
    {
        mutation_type = 1;
    }
    else
    {
        // Here we cannot remove a SNP
        if (pattern.size()==_pat_size_min)
        {
            uniform_int_distribution<int> distribution(0,1);
            // Choose the mutation type to perform (add/remove/change)
            mutation_type = distribution(_rng);
        }
        else
        {
            // Here we cannot add a SNP
            if (pattern.size()==_pat_size_max)
            {
                uniform_int_distribution<int> distribution(1,2);
                // Choose the mutation type to perform (add/remove/change)
                mutation_type = distribution(_rng);

            }
            // Here we can add, remove or change a SNP
            else
            {
                uniform_int_distribution<int> distribution(0,2);
                // Choose the mutation type to perform (add/remove/change)
                mutation_type = distribution(_rng);

            }
        }
    }

    switch (mutation_type)
    {
        // This case will add an snp to the pattern
        case 0:
        {
             // Pick a random SNP
            unsigned new_SNP = distribution_snp(_rng);
            // Repick while the snp is already in the pattern
            while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
            {
                new_SNP = distribution_snp(_rng);
            }
            // Add the SNP
            pattern.push_back(new_SNP);
        }
        break;
        // This case will change an snp from the pattern
        case 1:
        {
            // Pick a random SNP to change
            unsigned SNP_to_change = distribution_pattern(_rng);
            unsigned new_SNP = distribution_snp(_rng);
            // Repick while the snp is already in the pattern
            while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
            {
                new_SNP = distribution_snp(_rng);
            }
            pattern[SNP_to_change]= new_SNP;
        }
        break;
        // This case will remove an snp to the pattern
        case 2:
        {
            unsigned SNP_to_remove = distribution_pattern(_rng);
            pattern.erase(pattern.begin()+SNP_to_remove);
        }
        break;
        default:
        break;
    }
    if (k>1)
    {
        k = k-1;
        pattern = shake(pattern, k);
    }
    return pattern;
}

//==============================================================================
// vns : save_local_optimum
//==============================================================================
void vns::save_local_optimum(vector<unsigned> & x, vector<float> & x_score)
{
    sort(x.begin(), x.end());
    std::cout << " SNP in pattern:" << '\n';
    for (auto test : x)
    {
        std::cout << test << ' ';
    }
    std::cout << " Associated p-value: " << x_score[1] << '\n';
    // Search for the pattern in the result map
    auto current_opti = _optimum_set.find(x);

    if(_optimum_set.end() != current_opti)
    {
        // If the pattern is found, incrementation of the occurences number
        current_opti->second[0] += 1;
    }
    else
    {
        if (x_score[1] < _alpha)
        {
            // Add the solution
            _optimum_set[x] = {1, x_score[0], x_score[1], x_score[2]};
        }
    }
}

//==============================================================================
// vns : write_result_file
//==============================================================================

void vns::write_result_file()
{
    std::cout << "Time of execution: " << _duration << " seconds" << endl;
    // Create the output file name
    size_t firstindex = _filename.find_last_of("/");
    string filename_without_extension = _filename.substr(firstindex+1, 5000);
    size_t lastindex = filename_without_extension.find_last_of(".");
    filename_without_extension = filename_without_extension.substr(0, lastindex);

    std::cout << "### VNS has finished please see results in: " << '\n' << _output_directory + _output_prefix + filename_without_extension + "_result_vns.txt" << '\n';
    // Create the output file
    ofstream output_file(_output_directory + _output_prefix + filename_without_extension + "_result_vns.txt");

    output_file << "# Result from vns \n";
    output_file << "# Pattern || occurences || chi2-score || p-value || unreliable case\n";

    vector<pair<vector<unsigned>, vector<float>>> _optimum_set_vector;
    for (auto pattern : _optimum_set)
    {
        pair<vector<unsigned>, vector<float>> pair_sort(pattern.first, pattern.second);
        _optimum_set_vector.push_back(pair_sort);
    }
    sort(_optimum_set_vector.begin(), _optimum_set_vector.end(), compareFunc);

    for (auto const& pattern : _optimum_set_vector)
    {
        output_file << "{";
        for (auto const& snp : pattern.first)
        {
            output_file << _snp_id(snp);
            if (snp!=pattern.first.back())
            {
                output_file << ",";
            }
        }
        output_file << "}";
        unsigned cpt = 0;
        for (auto stat : pattern.second)
        {
            output_file << " || " << stat;
        }
        output_file << "\n";
    }
    output_file << "# Execution time : " << _duration << " seconds" << endl;
}
//==============================================================================
// vns : compareFunc
//==============================================================================
// Function to sort the results
bool vns::compareFunc(pair<vector<unsigned>, vector<float>> const& a, pair<vector<unsigned>, vector<float>> const& b)
{
    if (a.second[2] < b.second[2])
    {
        return true;
    }
    else
    {
        return false;
    }
}

//==============================================================================
// vns : print_parameters
//==============================================================================
void vns::print_parameters()
{
    std::cout << "### Parameters used for this run: " << endl;
    std::cout << "Number of iterations in VNS: " << _iteration_num  << endl;
    std::cout << "Alpha type I error rate: " << _alpha  << endl;
    std::cout << "Maximum size of pattern: " << _pat_size_max  << endl;
    std::cout << "Minimum size of pattern: " << _pat_size_min  << endl;
    std::cout << "k Max (vns stop counter) " << _k_max  << endl;
    std::cout << "l Max (local search stop counter) " << _l_max  << endl;
    std::cout << "### End of parameters" << endl;
    std::cout << endl;
}
