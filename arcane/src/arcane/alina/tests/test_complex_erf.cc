#define BOOST_TEST_MODULE TestComplex
#include <boost/test/unit_test.hpp>

#include <complex>

#include <amgcl/backend_builtin.h>
#include <amgcl/amg.h>
#include <amgcl/make_solver.h>

#include <amgcl/coarsening_aggregation.h>
#include <amgcl/coarsenin_smoothed_aggregation.h>
#include <amgcl/coarsening_smoothed_aggr_emin.h>

#include <amgcl/relaxation_damped_jacobi.h>
#include <amgcl/relaxation_gauss_seidel.h>
#include <amgcl/relaxation_spai0.h>
#include <amgcl/relaxation_ilu0.h>
#include <amgcl/relaxation_ilut.h>
#include <amgcl/relaxation_chebyshev.h>

#include <amgcl/solver_cg.h>
#include <amgcl/solver_bicgstab.h>
#include <amgcl/solver_gmres.h>

#include <amgcl/adapter_crs_tuple.h>
#include <amgcl/adapter_complex.h>
#include <amgcl/profiler.h>

#include "sample_problem.h"

namespace amgcl {
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

    typedef amgcl::backend::builtin<double> Backend;

    boost::property_tree::ptree prm;
    prm.put("precond.coarsening.aggr.block_size", 2);

    amgcl::make_solver<
        amgcl::amg<
            Backend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
        amgcl::solver::bicgstab<Backend>
        > solve( amgcl::adapter::complex_matrix(std::tie(n, ptr, col, val)), prm );

    std::cout << solve.precond() << std::endl;

    boost::iterator_range<const double*> f_range =
        amgcl::adapter::complex_range(rhs);

    boost::iterator_range<double*> x_range =
        amgcl::adapter::complex_range(x);

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
