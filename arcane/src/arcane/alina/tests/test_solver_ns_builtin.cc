#include <gtest/gtest.h>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/StaticMatrix.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_nonscalar_backend)
{
  test_backend< Alina::backend::BuiltinBackend< Alina::static_matrix<double, 2, 2> > >();
}
