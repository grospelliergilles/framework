#define BOOST_TEST_MODULE TestSolvers
#include <boost/test/unit_test.hpp>

#include <amgcl/value_type/complex.h>
#include <amgcl/backend_builtin.h>

#include "test_solver.h"

BOOST_AUTO_TEST_SUITE( test_solvers )

BOOST_AUTO_TEST_CASE(test_builtin_complex_backend)
{
    test_backend< amgcl::backend::builtin< std::complex<double> > >();
}

BOOST_AUTO_TEST_SUITE_END()
