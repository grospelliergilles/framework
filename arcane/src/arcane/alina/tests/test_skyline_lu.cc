#define BOOST_TEST_MODULE TestSkylineLU
#include <boost/test/unit_test.hpp>

#include <arcane/alina/adapter_zero_copy.h>
#include <arcane/alina/solver_skyline_lu.h>
#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/profiler.h>
#include "sample_problem.h"

namespace Arcane::Alina {
    profiler<> prof;
}
using namespace Arcane;

BOOST_AUTO_TEST_SUITE( test_skyline_lu )

BOOST_AUTO_TEST_CASE(skyline_lu)
{
    std::vector<ptrdiff_t> ptr;
    std::vector<ptrdiff_t> col;
    std::vector<double>    val;
    std::vector<double>    rhs;

    size_t n = sample_problem(16, val, col, ptr, rhs);

    auto A = Alina::adapter::zero_copy(n, ptr.data(), col.data(), val.data());

    Alina::solver::skyline_lu<double> solve(*A);

    std::vector<double> x(n);
    std::vector<double> r(n);

    solve(rhs, x);

    Alina::backend::residual(rhs, *A, x, r);

    BOOST_CHECK_SMALL(sqrt(Alina::backend::inner_product(r, r)), 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

