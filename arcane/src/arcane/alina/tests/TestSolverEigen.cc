// Pour Eigen
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include <gtest/gtest.h>

#include <arcane/alina/EigenBackend.h>

#include "TestSolverCommon.h"

TEST(alina_test_solvers, test_eigen_backend)
{
  test_backend<Alina::backend::EigenBackend<double>>();
}
