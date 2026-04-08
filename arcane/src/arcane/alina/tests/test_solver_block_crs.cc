#include <gtest/gtest.h>

#include <arcane/alina/BlockCSRBackend.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_block_crs_backend)
{
  test_backend< Alina::backend::BlockCSRBackend<double> >();
}
