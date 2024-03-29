/*
Authors: Tanguy Lallemand M2BB
         Jonathan Cruard M2BB
*/
#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random/discrete_distribution.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

//==============================================================================
// Creation of boost aliases
//==============================================================================
typedef boost::numeric::ublas::matrix<int> boost_matrix;
typedef boost::numeric::ublas::matrix<float> boost_matrix_float;
typedef boost::numeric::ublas::vector<float> boost_vector_float;
typedef boost::numeric::ublas::vector<int> boost_vector_int;
typedef boost::numeric::ublas::vector<std::string> boost_vector_string;
typedef boost::numeric::ublas::matrix_column<boost::numeric::ublas::matrix<int>> boost_column;
#endif
