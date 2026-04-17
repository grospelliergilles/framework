#include <gtest/gtest.h>

// To remove warnings about deprecated Eigen usage.
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include <Eigen/SparseLU>
#include "arccore/alina/EigenSolver.h"
#include "arccore/alina/BuiltinBackend.h"
#include "arccore/alina/Adapters.h"
#include "arccore/alina/Profiler.h"
#include "SampleProblemCommon.h"

namespace
{
Arcane::Alina::Profiler prof;
}

using namespace Arcane;

TEST(alina_test_solvers, eigen_solver)
{
  std::vector<int> ptr;
  std::vector<int> col;
  std::vector<double> val;
  std::vector<double> rhs;

  size_t n = sample_problem(16, val, col, ptr, rhs);
  Alina::CSRMatrix<double> A(std::tie(n, ptr, col, val));

  using Solver = Alina::EigenSolver<Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor, int>>>;

  Solver solve(A);

  std::vector<double> x(n);
  std::vector<double> r(n);

  solve(rhs, x);

  Alina::backend::residual(rhs, A, x, r);

  ASSERT_NEAR(sqrt(Alina::backend::inner_product(r, r)), 0.0, 1e-8);
}
