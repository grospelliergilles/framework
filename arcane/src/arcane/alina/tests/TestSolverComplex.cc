#include <gtest/gtest.h>

#include <arcane/alina/ValueTypeComplex.h>
#include <arcane/alina/BuiltinBackend.h>

#include "TestSolverCommon.h"

TEST(alina_test_solvers, test_builtin_complex_backend)
{
  test_backend< Alina::BuiltinBackend< std::complex<double> > >();
}
