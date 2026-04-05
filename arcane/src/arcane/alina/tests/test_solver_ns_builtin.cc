#define BOOST_TEST_MODULE TestSolvers
#include <boost/test/unit_test.hpp>
#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/value_type_static_matrix.h>

#include "test_solver.h"

BOOST_AUTO_TEST_SUITE( test_solvers )

BOOST_AUTO_TEST_CASE(test_nonscalar_backend)
{
    test_backend< amgcl::backend::builtin< amgcl::static_matrix<double, 2, 2> > >();
}

BOOST_AUTO_TEST_SUITE_END()
