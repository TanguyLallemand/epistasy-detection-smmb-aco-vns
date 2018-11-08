#include "contingency.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace std;

contingency::contingency(boost_matrix _genos_matrix, boost_matrix _phenos_matrix) /*: blas_dmatrix(2, 3)*/
{
    /*// Initialisation to zeros
    for (unsigned i=0; i < size1(); ++i)
    {
        for (unsigned j=0; j < size2(); ++j)
            this->at_element(i, j) = 0;
    }

    // Fill contingency table
    for(unsigned i=0; i<phenos.size(); ++i)
    {
        int row_contingency = phenos(i);
        int col_contingency = genos(i);
        if((row_contingency != 0 && row_contingency != 1) || (col_contingency != 0 && col_contingency != 1 && col_contingency != 2))
        {
            cout << "Index error while building contingency table." << endl;
            continue;
        }
        this->at_element(row_contingency, col_contingency) += 1;
    }*/
}
