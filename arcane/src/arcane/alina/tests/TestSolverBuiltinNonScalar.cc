#include <gtest/gtest.h>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/StaticMatrix.h>

#include "TestSolverCommon.h"

TEST(alina_test_solvers, test_nonscalar_backend)
{
  test_backend< Alina::BuiltinBackend< Alina::StaticMatrix<double, 2, 2> > >();
}
