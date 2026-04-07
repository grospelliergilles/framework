#include <gtest/gtest.h>

#include <arcane/alina/EigenBackend.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_eigen_backend)
{
  test_backend<Alina::backend::EigenBackend<double>>();
}
