#include <iostream>
#include <vector>

#include <../AlinaLib.h>
#include "sample_problem.h"

int main()
{
    std::vector<int>    ptr;
    std::vector<int>    col;
    std::vector<double> val;
    std::vector<double> rhs;

    int n = sample_problem(12l, val, col, ptr, rhs);

    amgclHandle prm = ARCANE_ALINA_params_create();

    ARCANE_ALINA_params_seti(prm, "precond.coarse_enough", 1000);
    ARCANE_ALINA_params_sets(prm, "precond.coarsening.type", "smoothed_aggregation");
    ARCANE_ALINA_params_setf(prm, "precond.coarsening.aggr.eps_strong", 1e-3f);
    ARCANE_ALINA_params_sets(prm, "precond.relax.type", "spai0");

    ARCANE_ALINA_params_sets(prm, "solver.type", "bicgstabl");
    ARCANE_ALINA_params_seti(prm, "solver.L", 1);
    ARCANE_ALINA_params_seti(prm, "solver.maxiter", 100);

    amgclHandle solver = ARCANE_ALINA_solver_create(
            n, ptr.data(), col.data(), val.data(), prm
            );

    ARCANE_ALINA_params_destroy(prm);

    std::vector<double> x(n, 0);
    conv_info cnv = ARCANE_ALINA_solver_solve(solver, rhs.data(), x.data());

    // Solve same problem again, but explicitly provide the matrix this time:
    std::fill(x.begin(), x.end(), 0);
    cnv = ARCANE_ALINA_solver_solve_mtx(
            solver, ptr.data(), col.data(), val.data(),
            rhs.data(), x.data()
            );

    std::cout << "Iterations: " << cnv.iterations << std::endl
              << "Error:      " << cnv.residual   << std::endl;

    ARCANE_ALINA_solver_destroy(solver);
}
