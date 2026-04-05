#define BOOST_TEST_MODULE TestComplex
#include <boost/test/unit_test.hpp>

#include <complex>

#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/amg.h>
#include <arcane/alina/make_solver.h>

#include <arcane/alina/coarsening_aggregation.h>
#include <arcane/alina/coarsening_smoothed_aggregation.h>
#include <arcane/alina/coarsening_smoothed_aggr_emin.h>

#include <arcane/alina/relaxation.h>

#include <arcane/alina/solver_cg.h>
#include <arcane/alina/solver_bicgstab.h>
#include <arcane/alina/solver_gmres.h>

#include <arcane/alina/adapter_crs_tuple.h>
#include <arcane/alina/adapter_complex.h>
#include <arcane/alina/profiler.h>

#include "sample_problem.h"

using namespace Arcane;

namespace Arcane::Alina {
    profiler<> prof;
}

BOOST_AUTO_TEST_SUITE( test_complex )

BOOST_AUTO_TEST_CASE(complex_matrix_adapter)
{
    typedef std::complex<double> complex;

    std::vector<int>     ptr;
    std::vector<int>     col;
    std::vector<complex> val;
    std::vector<complex> rhs;

    size_t n = sample_problem(32, val, col, ptr, rhs);

    std::vector<complex> x(n, complex(0.0,0.0));

    typedef Alina::backend::builtin<double> Backend;

    boost::property_tree::ptree prm;
    prm.put("precond.coarsening.aggr.block_size", 2);

    Alina::make_solver<
        Alina::amg<
            Backend,
            Alina::coarsening::smoothed_aggregation,
            Alina::relaxation::spai0
            >,
        Alina::solver::bicgstab<Backend>
        > solve( Alina::adapter::complex_matrix(std::tie(n, ptr, col, val)), prm );

    std::cout << solve.precond() << std::endl;

    boost::iterator_range<const double*> f_range =
        Alina::adapter::complex_range(rhs);

    boost::iterator_range<double*> x_range =
        Alina::adapter::complex_range(x);

    size_t iters;
    double resid;

    std::tie(iters, resid) = solve(f_range, x_range);

    BOOST_CHECK_SMALL(resid, 1e-8);

    std::cout
        << "iters: " << iters << std::endl
        << "resid: " << resid << std::endl
        ;
}

BOOST_AUTO_TEST_SUITE_END()
