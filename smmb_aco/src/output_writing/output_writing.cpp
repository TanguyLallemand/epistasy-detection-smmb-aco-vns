#include "output_writing.hpp"


void output_writing::write_in_file(string method_used)
{
    // Instanciate ofstream and open results.txt
    ofstream output_data("./results/results.txt");
    // // Open the file
    // output_data.open("./results/results.txt");
    // Check if file has trouble to be opened
    if( !output_data )
    {
      cerr << "Error: file could not be opened" << endl;
      // Return an error code
      exit(1);
    }
    // Output an header
    output_data << "Method used: " << method_used << endl;
    boost_vector_int num;
    for (size_t i = 0; i < 5/*taille amtrix*/; i++)
    {
        // Output datas
        output_data << num(i) << endl;
    }
    output_data.close();
}
