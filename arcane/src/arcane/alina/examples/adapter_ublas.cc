#include <iostream>
#include <cstdlib>
#include <utility>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include <arcane/alina/AMG.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/adapter_ublas.h>
#include <arcane/alina/coarsening.h>
#include <arcane/alina/relaxation.h>
#include <arcane/alina/solver_bicgstabl.h>
#include <arcane/alina/profiler.h>

#include "sample_problem.h"
using namespace Arcane;


typedef boost::numeric::ublas::compressed_matrix<
    double, boost::numeric::ublas::row_major
    > ublas_matrix;

typedef boost::numeric::ublas::vector<double> ublas_vector;

namespace Arcane::Alina {
    profiler<> prof;
}

int main(int argc, char *argv[]) {
    using Alina::prof;

    std::vector<int>    ptr;
    std::vector<int>    col;
    std::vector<double> val;
    std::vector<double> rhs;

    prof.tic("assemble");
    int m = argc > 1 ? atoi(argv[1]) : 128;
    int n = sample_problem(m, val, col, ptr, rhs);

    // Create ublas matrix with the data.
    ublas_matrix A(n, n);
    A.reserve(ptr[n]);

    for(int i = 0; i < n; ++i)
        for(int j = ptr[i], e = ptr[i+1]; j < e; ++j)
            A.push_back(i, col[j], val[j]);
    prof.toc("assemble");

    prof.tic("build");
    Alina::make_solver<
        Alina::AMG<
            Alina::backend::builtin<double>,
            Alina::coarsening::smoothed_aggregation,
            Alina::relaxation::spai0
            >,
        Alina::solver::BiCGStabLSolver<
            Alina::backend::builtin<double>
            >
        > solve( Alina::backend::map(A) );
    prof.toc("build");

    std::cout << solve.precond() << std::endl;

    ublas_vector x(n, 0);

    prof.tic("solve");
    Alina::SolverResult r = solve(rhs, x);
    prof.toc("solve");

    std::cout << "Iterations: " << r.nbIteration() << std::endl
              << "Error:      " << r.residual() << std::endl
              << std::endl << prof << std::endl;
}
