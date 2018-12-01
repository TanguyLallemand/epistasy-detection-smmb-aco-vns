#include "vns.hpp"

vns::vns(data_parsing dataset, parameters_parsing _params)
{
    //params unpacking
    this->_n_it_max = _params._n_it_max;
    this->_k_max = _params._k_max;
    this->_output_directory = _params.output_directory;
    this->_output_prefix = _params.output_prefix;

    //unpacking datas
    this->_genos_matrix = dataset._geno_matrix;
    this->_phenos_vector = dataset._pheno_vector;
    this->_snp_id = dataset._snp_id_vector;
    this->_filename = dataset._geno_filename.substr(14, dataset._geno_filename.length());

    this->_rng.seed(time(NULL));
}

//==============================================================================
// vns : run
//==============================================================================
void vns::run()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // Parallelization
    #pragma omp parallel for

    for (size_t i = 0; i < _n_it_max; i++)
    {
        std::cout << "iteration # : " << i << '\n';

        // Selecting starting pattern
        vector<unsigned> x = generate_starting_pattern();

        // compute the score of x for first iteration
        vector<float> x_score = test_pattern(x);

        // declaration of temporary stocking
        vector<unsigned> second_x;
        vector<unsigned> third_x;
        vector<float> third_x_score;

        // stop when x did not changed for k_max itérations
        int k = 0;
        while (k < _k_max)
        {
            // Take a random neighbor of x
            second_x = shake(x);

            // Searching for a best neighbor of second_x
            third_x_score = local_search(second_x, third_x);

            if (third_x_score[0] > x_score[0])
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
        // Saving the local optimum
        save_local_optimum(x, x_score);
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    double _duration = std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count();
    // Write results in a file
    write_result_file();
}


//==============================================================================
// vns : generate_starting_pattern
//==============================================================================
vector<unsigned> vns::generate_starting_pattern()
{
    // Initialization of list to store the pattern
    vector<unsigned> pattern;
    // Random pick of the pattern size up to _max_size (TODO fixed to 3 for now)
    std::uniform_int_distribution<int> distribution(1,3);
    unsigned size_pattern = distribution(_rng);

    // Fill the pattern with size_pattern different SNPs
    for (size_t i = 0; i < size_pattern; i++)
    {
        // Pick a random SNP
        std::uniform_int_distribution<int> distribution(0,_genos_matrix.size2()-1);
        unsigned new_SNP = distribution(_rng);

        while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
        {
            std::cout << "loop 1" << '\n';
            std::cout << new_SNP << '\n';
            for (auto tut : pattern)
            {
                std::cout << tut;
            }
            std::cout << '\n';
            new_SNP = distribution(_rng);
        }

        pattern.push_back(new_SNP);
    }

    return pattern;
}

//==============================================================================
// vns : local_search
//==============================================================================
vector<float> vns::local_search(vector<unsigned> second_x, vector<unsigned> & third_x)
{
    vector<float> score, best_score {0,0,0};
    unsigned k=0;
    vector<unsigned> candidat_neighbor;
    while (k<_k_max)
    {

        candidat_neighbor = shake(second_x);

        score = test_pattern(candidat_neighbor);

        if (score[0] > best_score[0])
        {
            best_score = score;
            third_x = candidat_neighbor;
            k=0;
        }
        else
        {
            k++;
        }
    }

    return best_score;
}

//==============================================================================
// vns : shake
//==============================================================================
vector<unsigned> vns::shake(vector<unsigned> pattern)
{
    unsigned mutation_type;
    std::uniform_int_distribution<int> distribution;
    std::uniform_int_distribution<int> distribution_pattern(0,pattern.size()-1);
    std::uniform_int_distribution<int> distribution_snp(0,_genos_matrix.size2()-1);
    // make sure we won't do a forbiden change in the pattern
    if (pattern.size()==1)
    {
        // Produce 0 or 1
        std::uniform_int_distribution<int> distribution(0,1);
    }
    else
    {
        // _max_size (TODO fixed to 3 for now)
        if (pattern.size()==3)
        {
            // Produce 1 or 2
            std::uniform_int_distribution<int> distribution(1,2);
        }
        else
        {
            // Produce 0, 1 or 2
            std::uniform_int_distribution<int> distribution(0,2);
        }
    }
    mutation_type = distribution(_rng);

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
                std::cout << "loop 2" << '\n';
                new_SNP = distribution_snp(_rng);
            }
            // add SNP
            pattern.push_back(new_SNP);
        }
        break;
        // This case will change an snp from the pattern
        case 1:
        {
            unsigned SNP_to_change = distribution_pattern(_rng);

            unsigned new_SNP = distribution_snp(_rng);
            // repick while the snp is already in the pattern
            while (find(pattern.begin(), pattern.end(), new_SNP) != pattern.end())
            {
                std::cout << "loop 3" << '\n';
                new_SNP = distribution_snp(_rng);
            }
            pattern[SNP_to_change]= new_SNP;
        }
        break;
        // This case will remove an snp to the pattern
        case 2: //delete
        {
            unsigned SNP_to_remove = distribution_pattern(_rng);

            pattern.erase(pattern.begin()+SNP_to_remove);
        }
        break;
        default :
        break;
    }

    return pattern;
}

//==============================================================================
// vns : save_local_optimum
//==============================================================================
void vns::save_local_optimum(vector<unsigned> & x, vector<float> & x_score)
{
    sort(x.begin(), x.end()); 
    auto current_opti = _optimum_set.find(x);
    if(_optimum_set.end() != current_opti)
    {
        current_opti->second[0] += 1;
    }
    else
    {
        _optimum_set[x] = {1, x_score[0], x_score[1], x_score[2]};
    }
}

//==============================================================================
// vns : save_local_optimum
//==============================================================================
void vns::write_result_file()
{
    std::cout << "Write_result_file" << '\n';
    std::cout << "Time of execution:" << _duration << "seconds" << endl;
    // Create the output file
    size_t lastindex = _filename.find_last_of(".");
    string filename_without_extension = _filename.substr(0, lastindex);
    std::cout << filename_without_extension << '\n';
    // Create the output file
    ofstream output_file(_output_directory + _output_prefix + filename_without_extension + "_vns.txt");

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
    output_file << "# Time of execution:" << _duration << "seconds" << endl;
}

//==============================================================================
// vns : test_pattern
//==============================================================================
vector<float> vns::test_pattern(vector<unsigned> const& pattern)
{
    vector<boost::numeric::ublas::matrix_column<boost_matrix>> pattern_datas;
    for (auto snp : pattern)
    {
        boost::numeric::ublas::matrix_column<boost_matrix> mc (_genos_matrix, snp);
        pattern_datas.push_back(mc);
    }
    vector<float> result = statistics::compute_p_value(pattern_datas, _phenos_vector);

    return result;
}
