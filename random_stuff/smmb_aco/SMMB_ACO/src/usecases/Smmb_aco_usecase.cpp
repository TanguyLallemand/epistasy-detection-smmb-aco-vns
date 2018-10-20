#include "Smmb_aco_usecase.hpp"
#include "Smmb_ACO.hpp"

#include <chrono>

using namespace std;
using namespace std::chrono;

Smmb_aco_usecase::Smmb_aco_usecase(blas_matrix & genos, blas_column & phenos, Parameters_file_parsing & params, blas_matrix & permuted_phenos)
    : _genos(genos), _phenos(phenos), _params(params)
{
    run(permuted_phenos);
}


void Smmb_aco_usecase::run(blas_matrix & permuted_phenos)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    Smmb_ACO smmb_ACO(_genos, _phenos, _params);
//    cout << "constructor of sm2b ok\n";
    smmb_ACO.run();
//    cout << "run of smb2b ok\n";
    smmb_ACO.make_consensus(permuted_phenos);
//    cout << "make_consensus of sm2b ok\n";

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(t2-t1).count();
    smmb_ACO.write_result(duration);
}
