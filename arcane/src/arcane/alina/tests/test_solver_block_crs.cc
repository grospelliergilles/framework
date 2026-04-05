#define BOOST_TEST_MODULE TestSolvers
#include <boost/test/unit_test.hpp>
#include <arcane/alina/backend_block_crs.h>

#include "test_solver.h"

BOOST_AUTO_TEST_SUITE( test_solvers )

BOOST_AUTO_TEST_CASE(test_block_crs_backend)
{
    test_backend< Alina::backend::block_crs<double> >();
}

BOOST_AUTO_TEST_SUITE_END()
