#include <vector>
#include <tuple>

#include <arcane/alina/adapter_crs_tuple.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/amg.h>
#include <arcane/alina/coarsenin_smoothed_aggregation.h>
#include <arcane/alina/relaxation_spai0.h>
#include <arcane/alina/solver_cg.h>
#include <arcane/alina/profiler.h>

#if defined(SOLVER_BACKEND_VEXCL)
#  include <amgcl/backend/vexcl.hpp>
   typedef amgcl::backend::vexcl<float>  fBackend;
   typedef amgcl::backend::vexcl<double> dBackend;
#else
#  ifndef SOLVER_BACKEND_BUILTIN
#    define SOLVER_BACKEND_BUILTIN
#  endif
#include <arcane/alina/backend_builtin.h>
typedef amgcl::backend::builtin<float>  fBackend;
   typedef amgcl::backend::builtin<double> dBackend;
#endif

#include "sample_problem.h"

   namespace amgcl { profiler<> prof; }
using amgcl::prof;

int main() {
    // Combine single-precision preconditioner with a
    // double-precision Krylov solver.
    typedef amgcl::make_solver<
        amgcl::amg<
            fBackend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
        amgcl::solver::cg< dBackend >
        >
        Solver;

    std::vector<ptrdiff_t> ptr, col;
    std::vector<double>    val, rhs;

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
    std::vector<double> &f = rhs;
    std::vector<double> x(n, 0.0);
#endif

    prof.tic("setup");
    Solver S(std::tie(n, ptr, col, val), Solver::params(), bprm);
    prof.toc("setup");

    std::cout << S << std::endl;

    int iters;
    double error;
    prof.tic("solve");
    std::tie(iters, error) = S(A_d, f, x);
    prof.toc("solve");

    std::cout << "Iterations: " << iters << std::endl
              << "Error:      " << error << std::endl
              << prof << std::endl;
}
