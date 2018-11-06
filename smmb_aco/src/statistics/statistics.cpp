#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <iostream>

#include "statistics.hpp"

//=================================================
// smmb_aco : constructeur //TODO work in progress
//=================================================
statistics::statistics(boost_matrix _genos_matrix, parameters_parsing _params)
{
    _alpha_phero = _params.aco_alpha;
    _beta_phero = _params.aco_beta;
    _tau = boost_vector(_genos_matrix.size2(), _params.aco_tau_init);
    _eta = boost_vector(_genos_matrix.size2(), _params.aco_eta);
}

//=================================================
// smmb_aco : generate_distribution
//=================================================
//Return un vecteur de distribution de probabilité
//TODO not tested
boost::numeric::ublas::vector<double> statistics::generate_distribution()
{
    boost::numeric::ublas::vector<double> probability;
    for (int i = 0; i < _tau.size(); i++) {
        probability.insert_element(probability.size(), _tau[i]*_alpha_phero + _eta[i]*_beta_phero);//si j ai bien compris c est l equivalent de pusb_back pour les vector de boost. rappel du proto insert_element(size_type  i, const_reference	t)
    };
    return probability;
}
