#include <gtest/gtest.h>

#include "arccore/alina/BlockCSRBackend.h"

#include "TestSolverCommon.h"

TEST(alina_test_solvers, test_block_crs_backend)
{
  test_backend< Alina::backend::BlockCSRBackend<double> >();
}
