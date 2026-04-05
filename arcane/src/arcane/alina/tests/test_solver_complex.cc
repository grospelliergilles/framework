#define BOOST_TEST_MODULE TestSolvers
#include <boost/test/unit_test.hpp>

#include <arcane/alina/value_type_complex.h>
#include <arcane/alina/backend_builtin.h>

#include "test_solver.h"

BOOST_AUTO_TEST_SUITE( test_solvers )

BOOST_AUTO_TEST_CASE(test_builtin_complex_backend)
{
    test_backend< Alina::backend::builtin< std::complex<double> > >();
}

BOOST_AUTO_TEST_SUITE_END()
