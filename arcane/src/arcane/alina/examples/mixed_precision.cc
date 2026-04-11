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
typedef Arcane::Alina::backend::BuiltinBackend<float> fBackend;
typedef Arcane::Alina::backend::BuiltinBackend<double> dBackend;

#include "sample_problem.h"

using namespace Arcane;
using namespace Arcane::Alina;

namespace
{
Profiler prof;
}

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

#ifdef SOLVER_BACKEND_VEXCL
  vex::Context ctx(vex::Filter::Env);
  std::cout << ctx << std::endl;

  bprm.q = ctx;
#endif

  prof.tic("assemble");
  int n = sample_problem(128, val, col, ptr, rhs);
  prof.toc("assemble");

#if defined(SOLVER_BACKEND_VEXCL)
  dBackend::matrix A_d(ctx, n, n, ptr, col, val);

  vex::vector<double> f(ctx, rhs);
  vex::vector<double> x(ctx, n);
  x = 0;
#elif defined(SOLVER_BACKEND_BUILTIN)
  auto A_d = std::tie(n, ptr, col, val);
  std::vector<double>& f = rhs;
  std::vector<double> x(n, 0.0);
#endif

  prof.tic("setup");
  Solver S(std::tie(n, ptr, col, val), Solver::params(), bprm);
  prof.toc("setup");

  std::cout << S << std::endl;

  prof.tic("solve");
  SolverResult r = S(A_d, f, x);
  prof.toc("solve");

  std::cout << "Iterations: " << r.nbIteration() << std::endl
            << "Error:      " << r.residual() << std::endl
            << prof << std::endl;
}
