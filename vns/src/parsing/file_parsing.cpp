/*
   Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
 */
 #include "file_parsing.hpp"

//=============================================================================
// data_parsing : constructor
//=============================================================================
data_parsing::data_parsing(string geno_filename, string pheno_filename, int header_size, char separator)
{
    _geno_filename = geno_filename;
    _pheno_filename = pheno_filename;
    _header_size = header_size;
    _separator = separator;
    get_line_nb();
    get_col_nb();

    //parsing each input object
    parse_snp_id();
    parse_geno();
    parse_pheno();
}

//=============================================================================
// data_parsing : parse_pheno
//=============================================================================
void data_parsing::parse_pheno()
{
    // Resize vector following data set size
    _pheno_vector.resize(_row_number);
    // Open file containing phenotype datas
    ifstream file(_pheno_filename);
    string line;
    for (size_t x = 0; x < _header_size; x++)
    {
        // Get header
        getline(file, line);
    }
    for (size_t i = 0; i < _pheno_vector.size(); i++)
    {
        // Get current line
        getline(file, line);
        // Save in right cell data, -48 to transform string to int using ASCII format
        _pheno_vector(i) = line[0]-48;
    }
}

//=============================================================================
// data_parsing : parse_geno
//=============================================================================
void data_parsing::parse_geno()
{
    // Rezise object to right dimensions
    _geno_matrix.resize( _row_number, _col_number);
    // Open file containing genos datas
    ifstream file(_geno_filename);
    // Init variable to handle line
    string line;
    // Get header data, storing SNPs names
    for (unsigned x = 0; x < _header_size; x++)
    {
        getline(file, line);
    }
    // Iterate tought file containing geno datas
    for (size_t i = 0; i < _geno_matrix.size1(); i++)
    {
        // Get line
        getline(file, line);
        // Init an iterator
        int r=0;
        // Iterate tought line
        for (size_t j = 0; j < line.size(); j++)
        {
            // Get only numbers
            if (line[j]!=',')
            {
                // Save in right cell data, -48 to transform string to int using ASCII format
                _geno_matrix(i, r) = line[j]-48;
                r++;
            }
        }
    }
}

//=================================================
// data_parsing : parse_snp_id
//=================================================
void data_parsing::parse_snp_id()
{
    _snp_id_vector.resize(_col_number);
    ifstream file(_geno_filename);
    string line;
    getline(file, line);
    int i=0;

    for (size_t r = 0; r < line.size(); r++)
    {
        if (line[r]==_separator)
        {
            i++;
        }
        else
        {
            _snp_id_vector(i) += line[r];
        }
    }
}
//=================================================
// data_parsing : get_line_nb
//=================================================
void data_parsing::get_line_nb()
{
    // Open file containing phenotype datas
    ifstream file(_pheno_filename);
    string temp;
    getline(file, temp);
    _row_number=0;
    // Check if file is opened
    if (file)
    {
        // Iterate tought line
        while (getline(file, temp))
        {
            // Increment row count
            _row_number++;
        }
    }
    else
    {
        cout << "unable to open : "<< _pheno_filename << '\n';
    }
}

//=================================================
// data_parsing : get_col_nb
//=================================================
void data_parsing::get_col_nb()
{
    // Open file containing genotype datas
    ifstream file(_geno_filename);
    string line;
    int i;
    _col_number=1;
    getline(file, line);
    int lenght = line.size();
    // Check if file is opened
    if (file)
    {
        // Iterate tought line
        for (i = 0; i < lenght; i++)
        {
            if (line[i]==_separator)
            {
                // Increment column count
                _col_number++;
            }
        }
    }
    else
    {
        cout << "unable to open : "<< _geno_filename << '\n';
    }
}
