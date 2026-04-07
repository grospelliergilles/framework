// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <vector>
#include <iostream>

#include <arcane/alina/backend_builtin_hybrid.h>
#include <arcane/alina/value_type_static_matrix.h>
#include <arcane/alina/Adapters.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/coarsening.h>
#include <arcane/alina/relaxation.h>
#include <arcane/alina/solver_cg.h>

#include <arcane/alina/IO.h>
#include <arcane/alina/profiler.h>

using namespace Arcane;

int main(int argc, char *argv[]) {
    // The command line should contain the matrix, the RHS, and the coordinate files:
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <A.mtx> <b.mtx> <coo.mtx>" << std::endl;
        return 1;
    }

    // The profiler:
    Alina::profiler<> prof("Nullspace");

    // Read the system matrix, the RHS, and the coordinates:
    ptrdiff_t rows, cols, ndim, ncoo;
    std::vector<ptrdiff_t> ptr, col;
    std::vector<double> val, rhs, coo;

    prof.tic("read");
    std::tie(rows, rows) = Alina::IO::mm_reader(argv[1])(ptr, col, val);
    std::tie(rows, cols) = Alina::IO::mm_reader(argv[2])(rhs);
    std::tie(ncoo, ndim) = Alina::IO::mm_reader(argv[3])(coo);
    prof.toc("read");

    Alina::precondition(ncoo * ndim == rows && (ndim == 2 || ndim == 3),
            "The coordinate file has wrong dimensions");

    std::cout << "Matrix " << argv[1] << ": " << rows << "x" << rows << std::endl;
    std::cout << "RHS "    << argv[2] << ": " << rows << "x" << cols << std::endl;
    std::cout << "Coords " << argv[3] << ": " << ncoo << "x" << ndim << std::endl;

    // Declare the solver type
    typedef Alina::static_matrix<double, 3, 3> DBlock;
    typedef Alina::static_matrix<float, 3, 3> FBlock;
    typedef Alina::backend::builtin_hybrid<DBlock> SBackend; // the solver backend
    typedef Alina::backend::builtin_hybrid<FBlock> PBackend; // the preconditioner backend

    typedef Alina::make_solver<
        Alina::AMG<
            PBackend,
            Alina::coarsening::smoothed_aggregation,
            Alina::relaxation::spai0
            >,
        Alina::solver::cg<SBackend>
        > Solver;

    // Solver parameters:
    Solver::params prm;
    prm.precond.coarsening.aggr.eps_strong = 0;

    // Convert the coordinates to the rigid body modes.
    // The function returns the number of near null-space vectors
    // (3 in 2D case, 6 in 3D case) and writes the vectors to the
    // std::vector<double> specified as the last argument:
    prm.precond.coarsening.nullspace.cols = Alina::coarsening::rigid_body_modes(
            ndim, coo, prm.precond.coarsening.nullspace.B);

    // We use the tuple of CRS arrays to represent the system matrix.
    auto A = std::tie(rows, ptr, col, val);

    // Initialize the solver with the system matrix.
    prof.tic("setup");
    Solver solve(A, prm);
    prof.toc("setup");

    // Show the mini-report on the constructed solver:
    std::cout << solve << std::endl;

    // Solve the system with the zero initial approximation:
    int iters;
    double error;
    std::vector<double> x(rows, 0.0);

    prof.tic("solve");
    std::tie(iters, error) = solve(A, rhs, x);
    prof.toc("solve");

    // Output the number of iterations, the relative error,
    // and the profiling data:
    std::cout << "Iters: " << iters << std::endl
              << "Error: " << error << std::endl
              << prof << std::endl;
}
