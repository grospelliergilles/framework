#include <gtest/gtest.h>

#include <complex>

#include "arccore/alina/BuiltinBackend.h"
#include "arccore/alina/AMG.h"
#include "arccore/alina/PreconditionedSolver.h"

#include "arccore/alina/Coarsening.h"
#include "arccore/alina/Relaxation.h"

#include "arccore/alina/ConjugateGradientSolver.h"
#include "arccore/alina/BiCGStabSolver.h"
#include "arccore/alina/GMRESSolver.h"

#include "arccore/alina/Adapters.h"
#include "arccore/alina/Profiler.h"

#include "SampleProblemCommon.h"

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

  typedef Alina::BuiltinBackend<double> Backend;

  Alina::PropertyTree prm;
  prm.put("precond.coarsening.aggr.block_size", 2);

  using AMGType = Alina::AMG<Backend, Alina::SmoothedAggregationCoarserning, Alina::SPAI0Relaxation>;
  using SolverType = Alina::PreconditionedSolver<AMGType, Alina::BiCGStabSolver<Backend>>;

  SolverType solve(Alina::adapter::complex_matrix(std::tie(n, ptr, col, val)), prm);

  std::cout << solve.precond() << std::endl;

  auto f_range = Alina::adapter::complex_range<const double>(rhs);
  auto x_range = Alina::adapter::complex_range<double>(x);

  SolverResult r = solve(f_range, x_range);

  ASSERT_NEAR(r.residual(), 0.0, 1e-8);

  std::cout << "iters: " << r.nbIteration() << std::endl
            << "resid: " << r.residual() << std::endl;
}
