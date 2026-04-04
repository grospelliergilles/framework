#define BOOST_TEST_MODULE TestSolvers
#include <boost/test/unit_test.hpp>
#include <amgcl/backend/EigenBackend.h>

#include "test_solver.h"

BOOST_AUTO_TEST_SUITE( test_solvers )

BOOST_AUTO_TEST_CASE(test_eigen_backend)
{
  test_backend<amgcl::backend::EigenBackend<double>>();
}

BOOST_AUTO_TEST_SUITE_END()
