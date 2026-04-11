#include <gtest/gtest.h>

#include <complex>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/PreconditionedSolver.h>

#include <arcane/alina/Coarsening.h>
#include <arcane/alina/Relaxation.h>

#include <arcane/alina/ConjugateGradientSolver.h>
#include <arcane/alina/BiCGStabSolver.h>
#include <arcane/alina/GMRESSolver.h>

#include <arcane/alina/Adapters.h>
#include <arcane/alina/Profiler.h>

#include "sample_problem.h"

using namespace Arcane;
using namespace Arcane::Alina;

namespace
{
Profiler prof;
}

TEST(alina_test_complex, complex_matrix_adapter)
{
  typedef std::complex<double> complex;

  std::vector<int> ptr;
  std::vector<int> col;
  std::vector<complex> val;
  std::vector<complex> rhs;

  size_t n = sample_problem(32, val, col, ptr, rhs);

  std::vector<complex> x(n, complex(0.0, 0.0));

  typedef Alina::backend::BuiltinBackend<double> Backend;

  Alina::PropertyTree prm;
  prm.put("precond.coarsening.aggr.block_size", 2);

  Alina::PreconditionedSolver<Alina::AMG<Backend,
                                Alina::SmoothedAggregationCoarserning,
                                Alina::SPAI0Relaxation>,
                     Alina::BiCGStabSolver<Backend>>
  solve(Alina::adapter::complex_matrix(std::tie(n, ptr, col, val)), prm);

  std::cout << solve.precond() << std::endl;

  boost::iterator_range<const double*> f_range = Alina::adapter::complex_range(rhs);
  boost::iterator_range<double*> x_range = Alina::adapter::complex_range(x);

  SolverResult r = solve(f_range, x_range);

  ASSERT_NEAR(r.residual(), 0.0, 1e-8);

  std::cout << "iters: " << r.nbIteration() << std::endl
            << "resid: " << r.residual() << std::endl;
}
