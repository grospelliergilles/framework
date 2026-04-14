// Pour Eigen
#pragma GCC diagnostic ignored "-Wdeprecated-copy"

#include <gtest/gtest.h>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/ValueTypeEigen.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_nonscalar_backend_eigen)
{
  test_backend< Alina::backend::BuiltinBackend< Eigen::Matrix<double, 2, 2> > >();
}
