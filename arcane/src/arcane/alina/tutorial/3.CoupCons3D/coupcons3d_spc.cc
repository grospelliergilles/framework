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

#include <iostream>
#include <string>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/Adapters.h>
#include <arcane/alina/value_type_static_matrix.h>
#include <arcane/alina/SchurPressureCorrectionPreconditioner.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/make_block_solver.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/BiCGStabSolver.h>
#include <arcane/alina/PreconditionerOnlySolver.h>
#include <arcane/alina/coarsening.h>
#include <arcane/alina/Relaxation.h>

#include <arcane/alina/IO.h>
#include <arcane/alina/profiler.h>

using namespace Arcane;

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    // The command line should contain the matrix file name:
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <matrix.mtx> <nu>" << std::endl;
        return 1;
    }

    // The profiler:
    Alina::profiler<> prof("CoupCons3D");

    // Read the system matrix:
    ptrdiff_t rows, cols;
    std::vector<ptrdiff_t> ptr, col;
    std::vector<double> val;

    prof.tic("read");
    std::tie(rows, cols) = Alina::IO::mm_reader(argv[1])(ptr, col, val);
    std::cout << "Matrix " << argv[1] << ": " << rows << "x" << cols << std::endl;
    prof.toc("read");

    // The RHS is filled with ones:
    std::vector<double> f(rows, 1.0);

    // The number of unknowns in the U subsystem
    ptrdiff_t nu = std::stoi(argv[2]);

    // We use the tuple of CRS arrays to represent the system matrix.
    // Note that std::tie creates a tuple of references, so no data is actually
    // copied here:
    auto A = std::tie(rows, ptr, col, val);

    // Compose the solver type
    typedef Alina::backend::BuiltinBackend<double> SBackend; // the outer iterative solver backend
    typedef Alina::backend::BuiltinBackend<float> PBackend;  // the PSolver backend
    typedef Alina::backend::BuiltinBackend<
        Alina::static_matrix<float,4,4>> UBackend;    // the USolver backend

    typedef Alina::make_solver<
        Alina::preconditioner::SchurPressureCorrectionPreconditioner<
            Alina::make_block_solver<
                Alina::AMG<
                    UBackend,
                    Alina::coarsening::aggregation,
                    Alina::relaxation::ILU0Relaxation
                    >,
                Alina::solver::PreconditionerOnlySolver<UBackend>
                >,
            Alina::make_solver<
                Alina::relaxation::as_preconditioner<
                    PBackend,
                    Alina::relaxation::SPAI0Relaxation
                    >,
                Alina::solver::PreconditionerOnlySolver<PBackend>
                >
            >,
        Alina::solver::BiCGStabSolver<SBackend>
        > Solver;

    // Solver parameters
    Solver::params prm;
    prm.precond.pmask.resize(rows);
    for(ptrdiff_t i = 0; i < rows; ++i) prm.precond.pmask[i] = (i >= nu);

    // Initialize the solver with the system matrix.
    prof.tic("setup");
    Solver solve(A, prm);
    prof.toc("setup");

    // Show the mini-report on the constructed solver:
    std::cout << solve << std::endl;

    // Solve the system with the zero initial approximation:
    std::vector<double> x(rows, 0.0);
    prof.tic("solve");
    Alina::SolverResult r = solve(A, f, x);
    prof.toc("solve");

    // Output the number of iterations, the relative error,
    // and the profiling data:
    std::cout << "Iters: " << r.nbIteration() << std::endl
              << "Error: " << r.residual() << std::endl
              << prof << std::endl;
}
