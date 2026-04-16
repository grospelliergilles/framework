#include <gtest/gtest.h>
#include <arcane/alina/BuiltinBackend.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_builtin_backend)
{
  test_backend< Alina::BuiltinBackend<double> >();
  test_backend< Alina::BuiltinBackend<double, int, ptrdiff_t> >();
  test_backend< Alina::BuiltinBackend<double, int, int> >();
  test_backend< Alina::BuiltinBackend<double, uint32_t, size_t> >();
}
