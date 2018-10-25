#include "Smmb_ACO.hpp"

#include <iostream>
#include <ctime>
#include <fstream>
#include <libgen.h>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "Services.hpp"
#include "G2_conditional_test_indep.hpp"
#include "Mutual_information.hpp"
#include "Permutations_adapt.hpp"

using namespace std;

//-----------------------------------------
// Smmb_ACO : Constructor
//-----------------------------------------
Smmb_ACO::Smmb_ACO(boost::numeric::ublas::matrix<int> & genos, 
					blas::matrix_column<blas::matrix<int> > & phenos, 
                    Parameters_file_parsing const& params)
					:  _genotypes(genos), _phenotypes(phenos), _params(params)
{
    // seed for randomness
    srand(time(NULL));
    _rng.seed(std::time(0));

    // Management of output files (traces, results).
    // File handlers will be open until SMMB destructor.
    string file_basename = basename((char*)_params.genos_file.c_str());
    string result_filename = "outputs/RESULT_" + file_basename;
    _results_handler.open(result_filename.c_str(), ios::trunc);

    if(!_results_handler)
    {
        std::cerr << "Error while opening output.txt (by writing access) !\n";
        exit(-1);
    }

    // Counter for the total number of tests of independance during a SMMB execution.
    _n_tests = 0;
    _tau_eta_sum = 0.0;
	
    unsigned n_snps = _genotypes.size2();

    _tau.assign(n_snps, (double)_params.aco_tau_init);
    _eta.assign(n_snps, (double)_params.aco_eta);
    _pdf.assign(n_snps, 0.0);
    _smmb_aco_exec_counter = params.n_smmb_aco_runs;
}

//-----------------------------------------
// Smmb : run
//-----------------------------------------
void Smmb_ACO::run()
{
    cout << "Running SMMB_ACO" << endl;

    // Learning multiple MBs with ACO

    for(unsigned aco_it_i=0 ; aco_it_i < _params.aco_n_iterations ; aco_it_i++)
    {
        cout << "ACO iteration " << aco_it_i << endl;
        tau_eta_sum();                  //cout << "Pheromone sum = " << _tau_eta_sum << endl;
        _scores_map.clear();

        // compute sampling proba for each snp, and corresponding cumulative distribution
        compute_aco_pdf();          //cout << "ACO PDF calculated\n";
        compute_aco_cumulative();   //cout << "ACO cumulative calculated\n";

        // Process all ants in the current iteration
        #pragma omp parallel for
        for(unsigned ant=0; ant<_params.aco_n_ants; ant++)
        {
            // select SNPs
            list<unsigned> snp_set;
            aco_sampling(snp_set);

            // learn MB from these SNPS
            list<unsigned> mb;
            learn_mb(mb, snp_set);

            // add candidate mb to _mbs
            if(!mb.empty())
                _mbs.push_back(mb);
        }

        // PDF is updated at each iteration in parallel synchronization
        update_tau();
    }
}

//-----------------------------------------------
// Smmb : run (after first SMMB-ACO execution)
//-----------------------------------------------
void Smmb_ACO::run(vector<unsigned> & snp_indexes)
{
    cout << "Running SMMB_ACO" << endl;

    // Learning multiple MBs with ACO
    // reinit PDF to match only snp_indexes
    // shorter vectors, corresponding to snps retained from previous iterations
    for(unsigned i=0; i<_tau.size(); i++)
    {
        if (find(snp_indexes.begin(), snp_indexes.end(), i) == snp_indexes.end())
            _tau[i] = 0;
    }

    for(unsigned aco_it_i=0 ; aco_it_i < _params.aco_n_iterations ; aco_it_i++)
    {
        cout << "ACO iteration " << aco_it_i << endl;
        tau_eta_sum();
        _scores_map.clear();

        // compute sampling proba for each snp, and corresponding cumulative distribution
        compute_aco_pdf();          //cout << "ACO PDF calculated\n";
        compute_aco_cumulative();   //cout << "ACO cumulative calculated\n";

        #pragma omp parallel for
        for(unsigned ant=0; ant<_params.aco_n_ants; ant++)
        {
            list<unsigned> snp_set;
            aco_sampling(snp_set);

            list<unsigned> mb;
            learn_mb(mb, snp_set);
            if(!mb.empty())
            {
//                cout << "mb added\n";
                _mbs.push_back(mb);
            }
        }
        update_tau();
    }
}

//-----------------------------------------
// Smmb_ACO : compute_aco_pdf
//-----------------------------------------
void Smmb_ACO::compute_aco_pdf()
{
    // P = f(pheromone)
    // vectors _tau and _eta have the same size (=nsnps)
    for(unsigned i=0; i<_tau.size(); i++)
        _pdf.at(i) = apply_aco_pdf_formula(_tau.at(i), _eta.at(i));
}

//-----------------------------------------
// Smmb_ACO : compute_aco_cumulative
//-----------------------------------------
void Smmb_ACO::compute_aco_cumulative()
{
    _cumulative_pdf.clear();
    double mem=0.0;
    unsigned i=0;
    for(auto& proba_i: _pdf)
    {
        if(proba_i > 0)
        {
            mem += proba_i;
            _cumulative_pdf[mem].push_back(i);
        }
        i++;
    }
}

//-----------------------------------------
// Smmb_ACO : aco_sampling
//-----------------------------------------
void Smmb_ACO::aco_sampling(list<unsigned> & snp_set)
{
    // to make sure that snp hasn't been already sampled
    bool first_time_sampled = false;

    unsigned set_size = _params.aco_set_size;
    for(unsigned i=0; i<set_size; i++)
    {
        // cout << "Filling " << i << " in set\n";
        // sampling number in [0;1[
        do
        {
            double sampled_proba = (double)rand()/(RAND_MAX);
            unsigned snp_index = get_snp_index_from_aco_cumulative(sampled_proba);

            if(find(snp_set.begin(), snp_set.end(), snp_index) == snp_set.end()) // not sampled yet
            {
                first_time_sampled = true;
                snp_set.push_back(snp_index);
            }
            else
            {
                first_time_sampled = false;
            }
        } while(! first_time_sampled);
    }
}

//-----------------------------------------
// Smmb_ACO : get_snp_index_from_aco_cumulative
//-----------------------------------------
unsigned Smmb_ACO::get_snp_index_from_aco_cumulative(double p)
{
//    cout << "Smmb_ACO::get_snp_index_from_aco_cumulative\n";
    auto it_up = _cumulative_pdf.lower_bound(p);
    auto it_low = it_up++;

    if(it_up != _cumulative_pdf.end())
    {
        auto it = it_low;
        if(it_low->first != p)
            it = it_up;

        if(it->second.size() == 1 )
        {
            // only 1 SNP with this proba to be sampled
            return it->second.front();
        }
        else
        {
            if(_smmb_aco_exec_counter < 2)
            {
                cout << "more than 2 snp with same key in cumulative pdf\n";
            }
            // more than 1 SNP with this proba to be sampled -> take 1 randomly
            return Services::draw_one(it->second, _rng);
        }
    }
    else
    {
        return 1;
    }
}

//-----------------------------------------
// Smmb_ACO : apply_aco_pdf_formula
//-----------------------------------------
double Smmb_ACO::apply_aco_pdf_formula(double pheromone_snp, double eta_snp)
{
    pheromone_snp = pow(pheromone_snp, _params.aco_alpha) * pow(eta_snp, _params.aco_beta);
    // sum of all probas = 1
    double proba_select = pheromone_snp / (double)_tau_eta_sum;

    return proba_select;
}

//-----------------------------------------
// Smmb_ACO : update_tau
//-----------------------------------------
void Smmb_ACO::update_tau()
{
    for(map<unsigned, list<double> >::const_iterator it=_scores_map.begin(); it!=_scores_map.end(); ++it)
    {
        unsigned snp_index = it->first;
        list<double> const& scores = it->second;

        for(list<double>::const_iterator list_it=scores.begin(); list_it!=scores.end(); list_it++)
        {
            update_tau(snp_index, *list_it);
        }
    }
}

//-----------------------------------------
// Smmb_ACO : update_tau
//-----------------------------------------
void Smmb_ACO::update_tau(unsigned snp_index, double assoc_score)
{
    _tau[snp_index] = (1-_params.aco_rho) * _tau[snp_index] + assoc_score*_params.aco_lambda;
}

//-----------------------------------------
// Smmb_ACO : tau_sum
//-----------------------------------------
void Smmb_ACO::tau_eta_sum()
{
    // Sum of pheromones over all SNPs
    _tau_eta_sum = 0.0;
    for(unsigned i=0; i<_tau.size(); i++)
        _tau_eta_sum += pow(_tau.at(i), _params.aco_alpha) * pow(_eta.at(i), _params.aco_beta);
}

//-----------------------------------------
// Smmb_ACO : learn_mb
//-----------------------------------------
void Smmb_ACO::learn_mb(list<unsigned> & mb, list<unsigned> & geno_indexes)
{
    unsigned cpt = 0;  // iteration counter for 1 MB learning
    list<unsigned> mem_mb; // useful to know if any changes happen during current iteration
    do  // learning of 1 MB
    {
        mem_mb = mb;
        forward(mb, geno_indexes);
        backward(mb, geno_indexes);
        cpt++;
    }while ((cpt < _params.n_trials_to_learn_1_mb && mb.empty()) && (mem_mb != mb));
    backward(mb, geno_indexes);
}

//-----------------------------------------
// Smmb_ACO : forward
//-----------------------------------------
void Smmb_ACO::forward(list<unsigned> & mb, list<unsigned> & geno_indexes)
{
    vector<unsigned> snps_drawnV;
    Services::random_subset(geno_indexes, snps_drawnV, _params.subset_size_small, _rng); // snps_drawn is given by access (reference)
    std::sort(snps_drawnV.begin(), snps_drawnV.end());                                                  // set sorted -> all subsets (powerset) will also be sorted.

    vector<vector<unsigned> > allsubsV;
    Services::powerset(snps_drawnV, allsubsV);

    // A vector containing all subsets (including empty set) of snps_drawn
//    vector<list<unsigned> > allsubs;
//    Services::powerset(snps_drawn, allsubs);

    // Analysis of the powerset of snps_drawn
    unsigned best_subset_index = 100;
    double  best_subset_pvalue = 1.1;

    for(unsigned i=1; i<allsubsV.size(); ++i) // Start at index 1 because of [0] --> empty list (due to powerset output)
    {
        vector<unsigned> const& current_combinV = allsubsV[i];   // reference to current subset
        list<unsigned> current_combin;                         // make a list from this vector
        for(unsigned const& i: current_combinV)
            current_combin.push_back(i);

//        list<unsigned> & current_combin = allsubs[i];   // reference to current subset

        /*  Explanation of below "for block"
         *  current_combination = {snp1, snp2, snp3}
         *  mb = {snp10}
         *  current_geno = snp1
         *  tmp_mb = {snp10, snp2, snp3}
         *  Indep_test (snp1, pheno | {snp10,snp2,snp3})    ===    Indep_test (current_geno, pheno | tmp_mb)
        */
        unsigned it_index = 0;
        for(list<unsigned>::const_iterator it = current_combin.begin(); it != current_combin.end(); it++)
        {
            unsigned current_geno = *it;
            list<unsigned> current_combin_trunc = current_combin;

            list<unsigned>::iterator iteratr = current_combin_trunc.begin();
            for(unsigned i=0; i<it_index; i++)
                iteratr ++;
            current_combin_trunc.erase(iteratr);

            // temporary merge mb and subset
            list<unsigned> tmp_mb = mb;
            Services::append_list_by_copy(tmp_mb, current_combin_trunc);

            // Conditional independance test
            G2_conditional_test_indep cond_g2(blas_column(_genotypes, current_geno), _phenotypes, tmp_mb);
            _n_tests ++;

            // A G2 independency test is reliable when there are enough observations in each cell of the contingency table.
            // See the method Contingency::is_reliable() for details.
            if(cond_g2.is_reliable())
            {
                if(cond_g2.pval() < best_subset_pvalue)
                {
                    best_subset_pvalue = cond_g2.pval();
                    best_subset_index = i;
                }

                #pragma omp critical
                _scores_map[i].push_back(cond_g2.g2());
            }
            it_index ++;
        }
    }
    // Append the best subset if alpha < threshold
    if(best_subset_pvalue < _params.alpha)
    {
        vector<unsigned> best_subset = allsubsV[best_subset_index];  //copy
        Services::append_vector_to_list(mb, best_subset);
        for(vector<unsigned>::iterator it=best_subset.begin(); it != best_subset.end(); ++it)
            geno_indexes.remove(*it);
    }
}

//-----------------------------------------
// Smmb_ACO : backward
//-----------------------------------------
void Smmb_ACO::backward(list<unsigned> & mb, list<unsigned> & geno_indexes)
{   
    list<unsigned>::iterator it = mb.begin();
    for(unsigned i=0; i<mb.size(); i++)
    {
        unsigned current = *it;
        vector<unsigned> mb_troncated;
        for(unsigned const& i: mb)
        {
            if(i != current)
                mb_troncated.push_back(i);
        }

        bool is_erased = false;
        vector<vector<unsigned> > allsubsets;
        Services::powerset(mb_troncated, allsubsets);
        for(unsigned j=0; j<allsubsets.size(); ++j)
        {
            list<unsigned> current_subset;
            for(unsigned const& elem: allsubsets[j])
                current_subset.push_back(elem);
            G2_conditional_test_indep cond_g2(blas_column(_genotypes, current), _phenotypes, current_subset);
            _n_tests++;

            if(cond_g2.pval() > _params.alpha) // non-reliable tests have assigned a p-value equal to 1
            {
                geno_indexes.push_front(*it);
                it = mb.erase(it);  // erase the current SNP from MB. Implicit iterator incrementation when erasing
                is_erased = true;
                break;              // if SNP is removed, we don't test it against other subsets of the Markov blanket
            }
        }
        if(!is_erased)
            ++it;
    }
}

//-----------------------------------------
// Smmb_ACO : make_consensus
//-----------------------------------------
void Smmb_ACO::make_consensus(blas_matrix & permuted_phenos)
{
    cout << "Running SMMB.make_consensus()" << endl;
    _consensus.clear();

    // Count occurrences of each Markov blanket and make them unique (mapping MB_list => number_of_occurences)
    map<list<unsigned>, unsigned> mb_occurrence;
    for(list<list<unsigned> >::iterator it=_mbs.begin(); it!=_mbs.end(); ++it)
    {
        list<unsigned>& current_mb = *it;
        current_mb.sort();
        auto mb_it = mb_occurrence.find(current_mb);
        if(mb_it == mb_occurrence.end())
            mb_occurrence[current_mb] = 1;
        else
            mb_occurrence[current_mb] += 1;
    }

    // Make a consensus (list object) of unique SNP ids
    for(map<list<unsigned>, unsigned>::const_iterator it = mb_occurrence.begin(); it != mb_occurrence.end(); ++it)
    {
        list<unsigned> const& li = it->first;
        Services::append_list_by_copy(_consensus, li);
    }

    // No duplicate SNP in the consensus
    _consensus.sort();
    _consensus.unique();

    if(_smmb_aco_exec_counter > 1 && _consensus.size() > _params.subset_size_small)
    {
        _smmb_aco_exec_counter--;
        _mbs.clear();

        vector<unsigned> consensus_snp_indexes;
        list<unsigned>::iterator it = _consensus.begin();
        for(unsigned i=0; i<_consensus.size(); i++)
        {
            consensus_snp_indexes.push_back(*it);
            it++;
        }

        // ajust aco_set_size (must not smaller than _consensus.size())
        if(_params.aco_set_size > _consensus.size())
            _params.aco_set_size = _consensus.size();

        run(consensus_snp_indexes);
        make_consensus(permuted_phenos);
    }
    else
    {
        // BACKWARD PHASE
        list<unsigned> mem_consensus;
        do
        {
            mem_consensus = _consensus;
            list<unsigned> removed_from_consensus;

            #pragma omp parallel for
            for(unsigned i=0; i<_consensus.size(); i++)
            {

                list<unsigned> updated_consensus =_consensus;
                Services::remove_listA_from_listB(removed_from_consensus, updated_consensus);

                list<unsigned>::iterator it = _consensus.begin();
                advance(it, i); // access with iterator to element to test for deletion
                unsigned current = *it;
                updated_consensus.remove(current);

                vector<unsigned> updated_consensusV;
                for(unsigned const& elem: updated_consensus)
                    updated_consensusV.push_back(elem);

                vector<vector<unsigned> > allsubsets;
                Services::powerset(updated_consensusV, allsubsets);

                for(unsigned j=0; j<allsubsets.size(); ++j)
                {
                    list<unsigned> current_subset;
                    for(unsigned const& elem: allsubsets[j])
                        current_subset.push_back(elem);

                    G2_conditional_test_indep cond_g2(blas_column(_genotypes, current), _phenotypes, current_subset);
                    if(cond_g2.is_reliable())
                    {
                        Permutations_adapt p_adapt(blas_column(_genotypes,current), _phenotypes, permuted_phenos, _params.alpha, _params.precision, "g2");
                        p_adapt.run();
                        cond_g2.set_pval(p_adapt.correction());
                    }

                    list<unsigned> keymap; keymap.push_back(*it);
                    Services::append_list_by_copy(keymap, current_subset);
                    _consensus_recorded_tests[keymap] = cond_g2.pval();

                    #pragma omp atomic
                    _n_tests ++;

                    if(cond_g2.pval() > _params.alpha)
                    {
                        #pragma omp critical
                        removed_from_consensus.push_front(current);
                        break;
                    }
                }
            }
            Services::remove_listA_from_listB(removed_from_consensus, _consensus);

        } while (_consensus != mem_consensus);
    }
}

//-----------------------------------------
// Smmb_ACO : write_result
//-----------------------------------------
void Smmb_ACO::write_result(double duration)
{
    _results_handler << "Duration\t" << duration << " milliseconds" << endl;
    _results_handler << "Output_file\t" << _params.genos_file << endl;
    _results_handler << "N_tests\t" << _n_tests << endl;
    _results_handler << "Markov_blanket consensus\t";
    if(!_consensus.empty())
    {
        list<unsigned>::const_iterator before_end=_consensus.end(); --before_end;
        _results_handler << "{ ";
        for(list<unsigned>::const_iterator it=_consensus.begin(); it!=before_end; ++it)
            _results_handler << *it << " ";
        _results_handler << _consensus.back() << " }" << endl;
    }
    else
    {
        _results_handler << "{ }";
    }
    _results_handler << endl << "### Occurences of weak learner MBs ####" << endl;
    map<list<unsigned>, unsigned> mb_occurrence;
    for(list<list<unsigned> >::iterator it=_mbs.begin(); it!=_mbs.end(); ++it)
    {
        list<unsigned>& current_mb = *it;
        current_mb.sort();
        auto mb_it = mb_occurrence.find(current_mb);
        if(mb_it == mb_occurrence.end())
            mb_occurrence[current_mb] = 1;
        else
            mb_occurrence[current_mb] += 1;
    }

    unsigned total_number_mbs = 0;
    for(map<list<unsigned>, unsigned>::const_iterator it = mb_occurrence.begin(); it != mb_occurrence.end(); ++it)
    {
        list<unsigned> const& li = it->first;
        _results_handler << "{ ";
        for (list<unsigned>::const_iterator l_it=li.begin(); l_it != li.end(); ++l_it)
            _results_handler << *l_it << " ";
        _results_handler << "} => " << it->second << "\n" ;
        total_number_mbs += it->second;
    }
    _results_handler << endl;
    _results_handler << "Total number of MBs learned : " << total_number_mbs << endl;

    _results_handler << endl << "### p_values ###" << endl;
    for(map<list<unsigned>, double>::const_iterator it = _consensus_recorded_tests.begin(); it != _consensus_recorded_tests.end(); ++it)
    {
        list<unsigned> const& current_combin = it->first;
        list<unsigned>::const_iterator l_it=current_combin.begin();

        _results_handler << "Indep( " << *l_it << ", phenotype ";
        if(current_combin.size() > 1)
        {
            _results_handler << "| ";
            advance(l_it, 1);
            for (; l_it != current_combin.end(); ++l_it)
                _results_handler << *l_it << " ";
            _results_handler << ") => " << it->second << "\n" ;
        }
        else
        {
            _results_handler << ")\n";
        }
    }
}
