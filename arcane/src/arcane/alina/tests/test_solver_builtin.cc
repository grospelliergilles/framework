#include <gtest/gtest.h>
#include <arcane/alina/backend_builtin.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_builtin_backend)
{
  test_backend< Alina::backend::builtin<double> >();
  test_backend< Alina::backend::builtin<double, int, ptrdiff_t> >();
  test_backend< Alina::backend::builtin<double, int, int> >();
  test_backend< Alina::backend::builtin<double, uint32_t, size_t> >();
}
