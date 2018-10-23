#include <list>
class smmb_aco
{
public:
void learn_MB(genotype_matrix, phenotype_matrix, int K, int n_it_n, float global_alpha, mem_a, P);
void forward(MB_modifie, MB_a, n_it_n, j, P, genotype_matrix, int K );
void backward(MB_a, phenotype_matrix, global_alpha);
private:
};
