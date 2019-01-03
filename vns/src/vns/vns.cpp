#include "vns.hpp"

//==============================================================================
// vns : vns
//==============================================================================
//Initialization of all variable for the run
vns::vns(data_parsing dataset, parameters_parsing _params)
{
//params unpacking
    //output informations
    this->_output_directory = _params.output_directory;
    this->_output_prefix = _params.output_prefix;
    //algorithm parameters
    this->_iteration_num = _params._iteration_num;
    this->_alpha = _params.alpha;
    this->_pat_size_max = _params._pat_size_max;
    this->_pat_size_min = _params._pat_size_min;
    this->_max_it_vns = _params._max_it_vns;
    this->_max_it_local_search = _params._max_it_local_search;
    //definition of neighbor exploration range
    this->_k_max = _params._pat_size_max-1;
    this->_l_max = _params._pat_size_max-1;

//unpacking datas
    this->_genos_matrix = dataset._geno_matrix;
    this->_phenos_vector = dataset._pheno_vector;
    this->_snp_id = dataset._snp_id_vector;
    this->_filename = dataset._geno_filename;

//Initialization of the rng seed for all random picks
    this->_rng.seed(time(NULL));
}

//==============================================================================
// vns : run
//==============================================================================
void vns::run()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    //Parallelization
    #pragma omp parallel for
    for (size_t i = 0; i < _iteration_num; i++)
    {
        std::cout << "iteration # : " << i << '\n';
        //Selecting starting pattern
        vector<unsigned> x = generate_starting_pattern();
        //compute the score of x for first iteration
        vector<float> x_score = test_pattern(x);
        //declaration of temporary stocking
        vector<unsigned> second_x, third_x;
        vector<float> third_x_score;

        //stop when x did not changed for k_max itérations
        unsigned exploration = 0;
        int k = 1;
        while (k < _k_max) //exploring neighborhood to a range of _k_max
        {
            while (exploration < _max_it_vns) //exploration of the current neighborhood
            {
                second_x = shake(x, k); //take a random neighbor of x
                third_x_score = local_search(second_x, third_x); //searching for a best neighbor of second_x
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
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    this-> _duration = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    std::cout << _duration << '\n';
    // Write results in a file
    write_result_file();
}


//==============================================================================
// vns : generate_starting_pattern
//==============================================================================
//take a pattern with a size between _pat_size_min and _pat_size_max and fill it
//with random SNPs
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
        //make sure to pick unique SNPs
        while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
        {
            new_SNP = distribution(_rng);
        }
        //add the selected SNP to the pattern
        pattern.push_back(new_SNP);
    }
    return pattern;
}

//==============================================================================
// vns : local_search
//==============================================================================
//exploration of the neighborhood of the provided pattern
vector<float> vns::local_search(vector<unsigned> second_x, vector<unsigned> & third_x)
{
    //Initialization of temporary variables
    vector<float> score, best_score {0,1,0};
    unsigned l=1;
    unsigned exploration = 0;
    while (l < _l_max)
    {
        while (exploration < _max_it_local_search)
        {
            vector<unsigned> candidat_neighbor = shake(second_x, l);
            score = test_pattern(candidat_neighbor);
            if (score[1] < best_score[1]) //if this pattern is the best tested he becomes the new best and reset counters
            {
                best_score = score;
                third_x = candidat_neighbor;
                l=1;
                exploration = 0;
            }
            else
            {
                exploration++;
            }
        }
        l++;
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
    //prepare the list of matrix column for the test
    vector<boost::numeric::ublas::matrix_column<boost_matrix>> pattern_datas;
    for (auto snp : pattern)
    {
        boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, snp);
        pattern_datas.push_back(mc);
    }
    //test the datas provided
    vector<float> result = statistics::compute_p_value(pattern_datas, _phenos_vector);
    return result;
}

//==============================================================================
// vns : shake
//==============================================================================
//return a neighbor of pattern, this neighbor have a distance of k with the pattern
vector<unsigned> vns::shake(vector<unsigned> pattern, unsigned k)
{
    unsigned mutation_type;
    std::uniform_int_distribution<int> distribution;
    std::uniform_int_distribution<int> distribution_pattern(0,pattern.size()-1);
    std::uniform_int_distribution<int> distribution_snp(0,_genos_matrix.size2()-1);
    // make sure we won't do a forbiden change in the pattern
    if (pattern.size()==_pat_size_min) //here we cannot remove a SNP
    {
        std::uniform_int_distribution<int> distribution(0,1);
    }
    else
    {
        if (pattern.size()==_pat_size_max) //here we cannot add a SNP
        {
            std::uniform_int_distribution<int> distribution(1,2);
        }
        else //here we can add, remove or change a SNP
        {
            std::uniform_int_distribution<int> distribution(0,2);
        }
    }
    //chose the mutation type to perform (add/remove/change)
    mutation_type = distribution(_rng);

    switch (mutation_type)
    {
        case 0: //this case will add an snp to the pattern
        {
            unsigned new_SNP = distribution_snp(_rng); // Pick a random SNP
            // Repick while the snp is already in the pattern
            while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
            {
                new_SNP = distribution_snp(_rng);
            }
            pattern.push_back(new_SNP); // add the SNP
        }
        break;
        case 1: //this case will change an snp from the pattern
        {
            unsigned SNP_to_change = distribution_pattern(_rng); //pick a random SNP to change
            unsigned new_SNP = distribution_snp(_rng);
            // repick while the snp is already in the pattern
            while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
            {
                new_SNP = distribution_snp(_rng);
            }
            pattern[SNP_to_change]= new_SNP;
        }
        break;
        case 2: //this case will remove an snp to the pattern
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
    //search for the pattern in the result map
    auto current_opti = _optimum_set.find(x);

    if(_optimum_set.end() != current_opti)
    {
        //if the pattern is found, incrementation of the occurences number
        current_opti->second[0] += 1;
    }
    else
    {
        if (x_score[1] < _alpha)
        {
            //add the solution
            _optimum_set[x] = {1, x_score[0], x_score[1], x_score[2]};
        }
    }
}

//==============================================================================
// vns : write_result_file
//==============================================================================
void vns::write_result_file()
{
    std::cout << "Write_result_file" << '\n';
    std::cout << "Time of execution: " << _duration << " seconds" << endl;
    // Create the output file name
    size_t firstindex = _filename.find_last_of("/");
    string filename_without_extension = _filename.substr(firstindex+1, 5000);
    size_t lastindex = filename_without_extension.find_last_of(".");
    filename_without_extension = filename_without_extension.substr(0, lastindex);

    std::cout << filename_without_extension << '\n';
    // Create the output file
    ofstream output_file(_output_directory + _output_prefix + filename_without_extension + "_result_vns.txt");

    output_file << "# Result from vns \n";
    output_file << "# Pattern || occurences || chi2-score || p-value || unreliable case\n";
    for (auto const& pattern : _optimum_set)
    {
        output_file << "{";
        for (auto const& snp : pattern.first)
        {
            std::cout << _snp_id(snp) << '\n';
            output_file << _snp_id(snp);
            if (snp!=pattern.first.back()) {
                output_file << ",";
            }
        }
        output_file << "}";
        for (auto stat : pattern.second)
        {
            output_file << " || " << stat;
        }
        output_file << "\n";
    }
    output_file << "# Execution time : " << _duration << " seconds" << endl;
}
