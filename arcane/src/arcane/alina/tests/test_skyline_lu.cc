#include <gtest/gtest.h>

#include <arcane/alina/Adapters.h>
#include <arcane/alina/SkylineLUSolver.h>
#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/Profiler.h>
#include "sample_problem.h"

namespace
{
Arcane::Alina::Profiler prof;
}

using namespace Arcane;

TEST(alina_test_skyline_lu, skyline_lu)
{
  std::vector<ptrdiff_t> ptr;
  std::vector<ptrdiff_t> col;
  std::vector<double> val;
  std::vector<double> rhs;

  size_t n = sample_problem(16, val, col, ptr, rhs);

  auto A = Alina::adapter::zero_copy(n, ptr.data(), col.data(), val.data());

  Alina::solver::SkylineLUSolver<double> solve(*A);

  std::vector<double> x(n);
  std::vector<double> r(n);

  solve(rhs, x);

  Alina::backend::residual(rhs, *A, x, r);

  ASSERT_NEAR(sqrt(Alina::backend::inner_product(r, r)), 0.0, 1e-8);
}
