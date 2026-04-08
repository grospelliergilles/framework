#include <gtest/gtest.h>

#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/value_type_eigen.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_nonscalar_backend_eigen)
{
  test_backend< Alina::backend::BuiltinBackend< Eigen::Matrix<double, 2, 2> > >();
}
