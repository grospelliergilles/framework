#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>

#include <arcane/alina/AMG.h>

#include <arcane/alina/Adapters.h>
#include <arcane/alina/backend_block_crs.h>
#include <arcane/alina/coarsening.h>
#include <arcane/alina/relaxation.h>
#include <arcane/alina/BiCGStabSolver.h>
#include <arcane/alina/profiler.h>

namespace Arcane::Alina
{
profiler<> prof("v2");
}
using namespace Arcane;

int main()
{
  using Alina::prof;

  typedef Alina::backend::BlockCSRBackend<double> Backend;
  typedef Alina::AMG<
  Backend,
  Alina::coarsening::aggregation,
  Alina::relaxation::spai0>
  AMG;

  std::vector<ptrdiff_t> ptr;
  std::vector<ptrdiff_t> col;
  std::vector<double> val;
  std::vector<double> rhs;

  prof.tic("read");
  {
    std::istream_iterator<int> iend;
    std::istream_iterator<double> dend;

    std::ifstream fptr("rows.txt");
    std::ifstream fcol("cols.txt");
    std::ifstream fval("values.txt");
    std::ifstream frhs("rhs.txt");

    Alina::precondition(fptr, "rows.txt not found");
    Alina::precondition(fcol, "cols.txt not found");
    Alina::precondition(fval, "values.txt not found");
    Alina::precondition(frhs, "rhs.txt not found");

    std::istream_iterator<int> iptr(fptr);
    std::istream_iterator<int> icol(fcol);
    std::istream_iterator<double> ival(fval);
    std::istream_iterator<double> irhs(frhs);

    ptr.assign(iptr, iend);
    col.assign(icol, iend);
    val.assign(ival, dend);
    rhs.assign(irhs, dend);
  }

  int n = ptr.size() - 1;
  prof.toc("read");

  prof.tic("build");
  AMG::params prm;
  prm.coarsening.aggr.eps_strong = 0;
  prm.coarsening.aggr.block_size = 4;
  prm.npre = prm.npost = 2;

  Backend::params bprm;
  bprm.block_size = 4;

  AMG amg(std::tie(n, ptr, col, val), prm, bprm);
  prof.toc("build");

  std::cout << amg << std::endl;

  std::vector<double> x(n, 0);

  Alina::solver::BiCGStabSolver<AMG::backend_type> solve(n);

  prof.tic("solve");
  Alina::SolverResult r = solve(amg, rhs, x);
  prof.toc("solve");

  std::cout << "Iterations: " << r.nbIteration() << std::endl
            << "Error:      " << r.residual() << std::endl;

  std::cout << Alina::prof << std::endl;
}
