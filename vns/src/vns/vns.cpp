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

    srand (time(NULL));
}

//==============================================================================
//vns : run
//==============================================================================
void vns::run()
{
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    // Initialisation of patterns and their neighbors
    generate_patterns();
    // Parallelization
    #pragma omp parallel for

    for (size_t i = 0; i < _n_it_max; i++)
    {
        std::cout << "iteration # : " << i << '\n';
        // Selecting starting pattern
        int index_pattern = rand() % (_pattern_list.size());

        // Initialization of starting pattern
        list<unsigned> x = _pattern_list[index_pattern];
        std::cout << "starting pattern" << '\n';
        for (auto test : x)
        {
            std::cout << test << ' ';
        }
        std::cout << '\n';
        // Initialization of vector to store g2 test's results
        vector<float> x_score(3);
        x_score = test_pattern(x);
        // Initialization of vector of list to store neighbors
        vector<list<unsigned>> x_neighbors;
        // Initialization of neighbors
        set_neighbors(x, x_neighbors);


        list<unsigned> second_x;
        vector<list<unsigned>> second_x_neighbors;

        list<unsigned> third_x;

        vector<float> third_x_score(3);

        int iterator = 0;

        int k = 0;
        while (k < _k_max)
        {
            // Take a random neighbor of x
            second_x = shake(x_neighbors);
            set_neighbors(second_x, second_x_neighbors);

            // Searching for the best neighbor of second_x
            third_x_score = variable_neighborhood_descent(second_x_neighbors, third_x);

            if (third_x_score[0] > x_score[0])
            {
                x = third_x;
                set_neighbors(x, x_neighbors);
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
    double _duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
    // Write results in a file
    write_result_file();
}


//==============================================================================
//vns : generate_patterns
//==============================================================================
void vns::generate_patterns()
{
    // Initialization of list to store SNP
    list<unsigned> temp;
    list<unsigned> snp_list(_genos_matrix.size2());
    iota(snp_list.begin(), snp_list.end(), 0);
    generate_patterns(temp, snp_list);
}

//==============================================================================
//vns : generate_patterns
//==============================================================================
void vns::generate_patterns(list<unsigned> temp, list<unsigned> snp_list)
{
    if (temp.size() < 3)
    {
        for (auto snp : snp_list)
        {
            // Add the snp to the pattern
            temp.push_back(snp);

            // Add the pattern to the vector
            _pattern_list.push_back(temp);

            // Copy the list of snp
            list<unsigned> next_snp_list = snp_list;

            // Remove snp added to temp pattern
            next_snp_list.remove(snp);

            // Recursive call to get next pattern
            generate_patterns(temp, next_snp_list);

            // Remove current snp from the list to let the next one come
            temp.pop_back();
        }
    }
    else
    {
        // Add the pattern to the vector
        _pattern_list.push_back(temp);
    }

}

//==============================================================================
//vns : set_neighbors
//==============================================================================
void vns::set_neighbors(list<unsigned> const& pattern, vector<list<unsigned>> & neighbors)
{
    neighbors.clear();
    //iterating all patterns to find neighbors for current at i+1 distance
    for (auto candidat_neighbor : _pattern_list)
    {
        //finding the number of common values needed to be neighbors at current distance
        unsigned neighbors_score = max(pattern.size(), candidat_neighbor.size())-1;

        //counting common values between vectors
        unsigned common_values = 0;
        for (auto s : pattern)
        {
            auto test = find(candidat_neighbor.begin(), candidat_neighbor.end(), s);
            if (*test != s)
            {
                common_values+=1;
            }
        }

        //if the 2 sets have only 1 difference they are neighbors
        if (neighbors_score == common_values)
        {
            //append the new neighbor to the list of neighbors with i+1 distance
            neighbors.push_back(candidat_neighbor);
        }
    }
}

//==============================================================================
//vns : variable_neighborhood_descent
//==============================================================================
vector<float> vns::variable_neighborhood_descent(vector<list<unsigned>> const& neighbors, list<unsigned> & third_x)
{
    vector<float> score, best_score = {0,0,0};

    for (auto neighbor_iterator : neighbors)
    {
        score = test_pattern(neighbor_iterator);

        if (score[0] > best_score[0])
        {
            best_score = score;
            third_x = neighbor_iterator;
        }
    }
    return best_score;
}

//==============================================================================
//vns : shake
//==============================================================================
list<unsigned> vns::shake(vector<list<unsigned>> neighbors)
{
    int index_pattern = rand() % (neighbors.size());

    list<unsigned> x2 = neighbors[index_pattern];

    return x2;
}

//==============================================================================
//vns : save_local_optimum
//==============================================================================
void vns::save_local_optimum(list<unsigned> & x, vector<float> & x_score)
{
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
//vns : save_local_optimum
//==============================================================================
void vns::write_result_file()
{
    std::cout << "Write_result_file" << '\n';
    std::cout << "Time of execution:" << _duration;
    //create the output file
    size_t lastindex = _filename.find_last_of(".");
    string filename_without_extension = _filename.substr(0, lastindex);
    std::cout << filename_without_extension << '\n';
    //create the output file
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
    output_file << "# Time of execution:" << _duration;
}

//==============================================================================
//vns : test_pattern
//==============================================================================
vector<float> vns::test_pattern(list<unsigned> const& pattern)
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
