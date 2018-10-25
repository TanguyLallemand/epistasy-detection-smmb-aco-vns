#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP

#include <list>
#include "Parameters_file_parsing.hpp"
class smmb_aco
{
public:

void learn_MB(genotype_matrix, phenotype_matrix, int K, int n_it_n, double alpha, int mem_a, P);
void forward(MB_modifie, MB_a, int n_it_n, int j, P, genotype_matrix, int K );
void backward(MB_a, phenotype_matrix, double alpha);
void run();

private:

int mem_a;
P;
double alpha;
MB_A;
int j;
int K;
bool MB_modified
S;
parameters_parsing _params;
};
#endif
