#include <vector>
#include <tuple>

#include <arcane/alina/Adapters.h>
#include <arcane/alina/PreconditionedSolver.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/Coarsening.h>
#include <arcane/alina/Relaxation.h>
#include <arcane/alina/ConjugateGradientSolver.h>
#include <arcane/alina/Profiler.h>

#ifndef SOLVER_BACKEND_BUILTIN
#define SOLVER_BACKEND_BUILTIN
#endif
#include <arcane/alina/BuiltinBackend.h>
typedef Arcane::Alina::BuiltinBackend<float> fBackend;
typedef Arcane::Alina::BuiltinBackend<double> dBackend;

#include "sample_problem.h"

using namespace Arcane;
using namespace Arcane::Alina;

int main()
{
  // Combine single-precision preconditioner with a
  // double-precision Krylov solver.
  typedef Alina::PreconditionedSolver<Alina::AMG<fBackend,
                                        Alina::SmoothedAggregationCoarserning,
                                        Alina::SPAI0Relaxation>,
                             Alina::ConjugateGradientSolver<dBackend>>
  Solver;

  std::vector<ptrdiff_t> ptr, col;
  std::vector<double> val, rhs;

  dBackend::params bprm;

  ARCANE_ALINA_TIC("assemble");
  int n = sample_problem(128, val, col, ptr, rhs);
  ARCANE_ALINA_TOC("assemble");

  auto A_d = std::tie(n, ptr, col, val);
  std::vector<double>& f = rhs;
  std::vector<double> x(n, 0.0);

  ARCANE_ALINA_TIC("setup");
  Solver S(std::tie(n, ptr, col, val), Solver::params(), bprm);
  ARCANE_ALINA_TIC("setup");

  std::cout << S << std::endl;

  ARCANE_ALINA_TIC("solve");
  SolverResult r = S(A_d, f, x);
  ARCANE_ALINA_TIC("solve");

  std::cout << "Iterations: " << r.nbIteration() << std::endl
            << "Error:      " << r.residual() << std::endl
            << Alina::Profiler::globalProfiler() << std::endl;
}
