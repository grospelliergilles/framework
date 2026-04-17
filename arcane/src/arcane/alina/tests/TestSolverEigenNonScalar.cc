// Pour Eigen
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include <gtest/gtest.h>

#include "arcane/alina/BuiltinBackend.h"
#include "arcane/alina/ValueTypeEigen.h"

#include "TestSolverCommon.h"

TEST(alina_test_solvers, test_nonscalar_backend_eigen)
{
  test_backend< Alina::BuiltinBackend< Eigen::Matrix<double, 2, 2> > >();
}
