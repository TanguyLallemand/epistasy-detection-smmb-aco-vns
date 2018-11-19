#include <iostream>
#include <fstream>
// For exit function
#include <cstdlib>

#include "global.hpp"
#include "output_writing.hpp"

using std::ofstream;
using std::cerr;
using std::endl;

void output_writing::write_in_file()
{
    ofstream output_data;
    // Open the file
    output_data.open("./results/results.txt");
    // Check if file has trouble to be opened
    if( !output_data )
    {
      cerr << "Error: file could not be opened" << endl;
      // Return an error code
      exit(1);
    }
    boost_vector_int num;
    for (size_t i = 0; i < 5/*taille amtrix*/; i++)
    {
        // Output datas
        output_data << num(i) << endl;
    }
    output_data.close();
}
