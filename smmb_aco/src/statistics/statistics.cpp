#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"

//=================================================
// smmb_aco : constructeur //TODO work in progress
//=================================================
statistics::statistics(boost_matrix _genos_matrix, boost_matrix _phenos_matrix, parameters_parsing _params)
{
    _alpha_phero = _params.aco_alpha;
    _beta_phero = _params.aco_beta;
    _tau = boost_vector(_genos_matrix.size2(), _params.aco_tau_init);
    _eta = boost_vector(_genos_matrix.size2(), _params.aco_eta);
}

//-----------------------------------------
// chi2_test_indep
//-----------------------------------------
statistics::chi2_test_indep()
{
    //need contignecy
    run();
}

//-----------------------------------------
// run a chi2_test_indep
//-----------------------------------------
statistics::run()
{

}
