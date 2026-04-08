#include <gtest/gtest.h>

#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/value_type_static_matrix.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_nonscalar_backend)
{
  test_backend< Alina::backend::BuiltinBackend< Alina::static_matrix<double, 2, 2> > >();
}
