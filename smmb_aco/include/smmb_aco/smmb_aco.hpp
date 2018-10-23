#include <list>
class smmb_aco
{
public:
void learn_MB(genotype_matrix, phenotype_matrix, int K, int n_it_n, float global_alpha, mem_a, P);
void forward(std::list<unsigned> & mb, std::list<unsigned> & genotype_indexes);
void backward(std::list<unsigned> & mb, std::list<unsigned> & genotype_indexes);
private:
};
