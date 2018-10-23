#include <list>
class smmb_aco
{
public:
void learn_MB(std::list<unsigned> & mb, std::list<unsigned> & genotype_indexes);
void forward(std::list<unsigned> & mb, std::list<unsigned> & genotype_indexes);
void backward(std::list<unsigned> & mb, std::list<unsigned> & genotype_indexes);
private:
};
