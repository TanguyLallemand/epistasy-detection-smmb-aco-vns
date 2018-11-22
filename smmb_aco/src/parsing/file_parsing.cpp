#include "file_parsing.hpp"

//TODO voir si un destructeur est nécessaire

//=============================================================================
// data_parsing : constructeur
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
// data_parsing : parse_geno
//=============================================================================
void data_parsing::parse_geno()
{
    _geno_matrix.resize( _row_number, _col_number);
    ifstream file(_geno_filename);
    string line;
    for (unsigned x = 0; x < _header_size; x++)
    {
        getline(file, line);
    }
    for (size_t i = 0; i < _geno_matrix.size1(); i++)
    {
        getline(file, line);
        int r=0;
        for (size_t j = 0; j < line.size(); j++)
        {
            if (line[j]!=',')
            {
                _geno_matrix(i, r) = line[j]-48;
                r++;
            }
        }
    }
}

//=============================================================================
// data_parsing : parse_geno
//=============================================================================
void data_parsing::parse_pheno()
{
    _pheno_vector.resize(_row_number);
    ifstream file(_pheno_filename);
    string line;
    for (size_t x = 0; x < _header_size; x++)
    {
        getline(file, line);
    }
    for (size_t i = 0; i < _pheno_vector.size(); i++)
    {
        getline(file, line);
        _pheno_vector(i) = line[0]-48;
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
    ifstream file(_pheno_filename);
    string temp;
    getline(file, temp);
    _row_number=0;
    if (file)
    {
        while (getline(file, temp))
        {
            _row_number++;
        }
        //_row_number = _row_number - _header_size;
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
    ifstream file(_geno_filename);
    string line;
    int i;
    _col_number=1;
    getline(file, line);
    int lenght = line.size();
    for (i = 0; i < lenght; i++)
    {
        if (line[i]==_separator)
        {
            _col_number++;
        }
    }
}
