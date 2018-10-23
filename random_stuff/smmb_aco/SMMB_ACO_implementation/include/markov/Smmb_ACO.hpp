#ifndef SMMB_ACO_HPP
#define SMMB_ACO_HPP

#include <list>
#include <map>
#include <utility>
#include <string>
#include <random>
#include <ctime>
#include <fstream>
#include <boost/numeric/ublas/matrix.hpp>

#include "common.h"
#include "Parameters_file_parsing.hpp"

class Smmb_ACO 
{
public:
	Smmb_ACO(boost::numeric::ublas::matrix<int> & genos, 
			 blas::matrix_column<blas::matrix<int> > & phenos, 
             Parameters_file_parsing const& params);

    void run();
    void run(std::vector<unsigned> & snp_indexes);

    unsigned get_snp_index_from_aco_cumulative(double p);
    void compute_aco_pdf();
    void compute_aco_cumulative();
    void aco_sampling(std::list<unsigned> & snp_set);

    double apply_aco_pdf_formula(double pheromone_snp, double eta_snp);
    void update_tau();
    void update_tau(unsigned snp_index, double assoc_score);

    void tau_eta_sum();

    void learn_mb(std::list<unsigned> & mb, std::list<unsigned> & geno_indexes);
    void forward(std::list<unsigned> & mb, std::list<unsigned> & geno_indexes);
    void backward(std::list<unsigned> & mb, std::list<unsigned> & geno_indexes);

    void make_consensus(blas_matrix & permuted_phenos);
    void write_result(double duration);
	
	
protected:
	boost::numeric::ublas::matrix<int> & _genotypes;
    blas::matrix_column<blas::matrix<int> > & _phenotypes;
    Parameters_file_parsing _params;
	std::list<std::list<unsigned> > _mbs;       // list of all learned MBs
    std::list<unsigned> _consensus;             // consensus of all learned MBs
    std::mt19937 _rng;                          // randomness generator
	std::ofstream _results_handler;             // output file handler
    std::ofstream _trace_handler;               // trace file handler
    unsigned _n_tests;
    double _tau_eta_sum;
    std::vector<double> _tau;					// probability distribution function
    std::vector<double> _eta;
    std::vector<double> _pdf;
    std::map<double, std::vector<unsigned> > _cumulative_pdf;     // cumulative probability distribution function
    std::map<unsigned, std::list<double> > _scores_map;
    std::map<std::list<unsigned>, double>  _consensus_recorded_tests;
    unsigned _smmb_aco_exec_counter;                // count number of executions of smmb
};

#endif
