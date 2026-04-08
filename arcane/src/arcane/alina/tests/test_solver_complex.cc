#include <gtest/gtest.h>

#include <arcane/alina/value_type_complex.h>
#include <arcane/alina/BuiltinBackend.h>

#include "test_solver.h"

TEST(alina_test_solvers, test_builtin_complex_backend)
{
  test_backend< Alina::backend::BuiltinBackend< std::complex<double> > >();
}
