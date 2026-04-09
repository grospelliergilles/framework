#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/value_type_static_matrix.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/solver_runtime.h>
#include <arcane/alina/CoarseningRuntime.h>
#include <arcane/alina/RelaxationRuntime.h>
#include <arcane/alina/relaxation.h>
#include <arcane/alina/CPRPreconditioner.h>
#include <arcane/alina/Adapters.h>
#include <arcane/alina/IO.h>
#include <arcane/alina/profiler.h>

using namespace Arcane;

namespace Arcane::Alina { profiler<> prof; }
using Alina::prof;
using Alina::precondition;

//---------------------------------------------------------------------------
template <class Matrix>
void solve_cpr(const Matrix& K, const std::vector<double>& rhs, Alina::PropertyTree& prm)
{
  auto t1 = prof.scoped_tic("CPR");

  using Backend = Alina::backend::BuiltinBackend<double>;
  using PPrecond = Alina::AMG<Backend, Alina::runtime::coarsening::CoarseningRuntime, Alina::runtime::relaxation::RuntimeRelaxation>;
  using SPrecond = Alina::relaxation::as_preconditioner<Backend, Alina::runtime::relaxation::RuntimeRelaxation>;

  prof.tic("setup");
  Alina::make_solver<Alina::preconditioner::CPRPreconditioner<PPrecond, SPrecond>,
                     Alina::runtime::solver::wrapper<Backend>>
  solve(K, prm);
  prof.toc("setup");

  std::cout << solve.precond() << std::endl;

  std::vector<double> x(rhs.size(), 0.0);

  prof.tic("setup");
  Alina::SolverResult r = solve(rhs, x);
  prof.toc("setup");

  std::cout << "Iterations: " << r.nbIteration() << std::endl
            << "Error:      " << r.residual() << std::endl;
}

//---------------------------------------------------------------------------
template <int B, class Matrix>
void solve_block_cpr(const Matrix& K, const std::vector<double>& rhs, Alina::PropertyTree& prm)
{
  auto t1 = prof.scoped_tic("CPR");

  typedef Alina::static_matrix<double, B, B> val_type;
  typedef Alina::static_matrix<double, B, 1> rhs_type;
  typedef Alina::backend::BuiltinBackend<val_type> SBackend;
  typedef Alina::backend::BuiltinBackend<double> PBackend;

  typedef Alina::AMG<
  PBackend,
  Alina::runtime::coarsening::CoarseningRuntime,
  Alina::runtime::relaxation::RuntimeRelaxation>
  PPrecond;

  typedef Alina::relaxation::as_preconditioner<
  SBackend,
  Alina::runtime::relaxation::RuntimeRelaxation>
  SPrecond;

  prof.tic("setup");
  Alina::make_solver<
  Alina::preconditioner::CPRPreconditioner<PPrecond, SPrecond>,
  Alina::runtime::solver::wrapper<SBackend>>
  solve(Alina::adapter::block_matrix<val_type>(K), prm);
  prof.toc("setup");

  std::cout << solve.precond() << std::endl;

  std::vector<rhs_type> x(rhs.size(), Alina::math::zero<rhs_type>());

  auto rhs_ptr = reinterpret_cast<const rhs_type*>(rhs.data());
  size_t n = Alina::backend::rows(K) / B;

  prof.tic("solve");
  Alina::SolverResult r = solve(Alina::make_iterator_range(rhs_ptr, rhs_ptr + n), x);
  prof.toc("solve");

  std::cout << "Iterations: " << r.nbIteration() << std::endl
            << "Error:      " << r.residual() << std::endl;
}

//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  using Alina::precondition;
  using Alina::prof;
  using std::string;
  using std::vector;

  namespace po = boost::program_options;
  namespace io = Alina::IO;

  po::options_description desc("Options");

  desc.add_options()("help,h", "show help")(
  "binary,B",
  po::bool_switch()->default_value(false),
  "When specified, treat input files as binary instead of as MatrixMarket. "
  "It is assumed the files were converted to binary format with mm2bin utility. ")(
  "matrix,A",
  po::value<string>()->required(),
  "The system matrix in MatrixMarket format")(
  "rhs,f",
  po::value<string>(),
  "The right-hand side in MatrixMarket format")(
  "runtime-block-size,b",
  po::value<int>(),
  "The block size of the system matrix set at runtime")(
  "static-block-size,c",
  po::value<int>()->default_value(1),
  "The block size of the system matrix set at compiletime")(
  "params,P",
  po::value<string>(),
  "parameter file in json format")(
  "prm,p",
  po::value<vector<string>>()->multitoken(),
  "Parameters specified as name=value pairs. "
  "May be provided multiple times. Examples:\n"
  "  -p solver.tol=1e-3\n"
  "  -p precond.coarse_enough=300");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  po::notify(vm);

  Alina::PropertyTree prm;
  if (vm.count("params"))
    prm.read_json(vm["params"].as<string>());

  if (vm.count("prm")) {
    for (const string& v : vm["prm"].as<vector<string>>()) {
      Alina::put(prm, v);
    }
  }

  int cb = vm["static-block-size"].as<int>();

  if (vm.count("runtime-block-size"))
    prm.put("precond.block_size", vm["runtime-block-size"].as<int>());
  else
    prm.put("precond.block_size", cb);

  size_t rows;
  vector<ptrdiff_t> ptr, col;
  vector<double> val, rhs;
  std::vector<char> pm;

  {
    auto t = prof.scoped_tic("reading");

    string Afile = vm["matrix"].as<string>();
    bool binary = vm["binary"].as<bool>();

    if (binary) {
      io::read_crs(Afile, rows, ptr, col, val);
    }
    else {
      size_t cols;
      std::tie(rows, cols) = io::mm_reader(Afile)(ptr, col, val);
      precondition(rows == cols, "Non-square system matrix");
    }

    if (vm.count("rhs")) {
      string bfile = vm["rhs"].as<string>();

      size_t n, m;

      if (binary) {
        io::read_dense(bfile, n, m, rhs);
      }
      else {
        std::tie(n, m) = io::mm_reader(bfile)(rhs);
      }

      precondition(n == rows && m == 1, "The RHS vector has wrong size");
    }
    else {
      rhs.resize(rows, 1.0);
    }
  }

#define CALL_BLOCK_SOLVER(z, data, B) \
  case B: \
    solve_block_cpr<B>(std::tie(rows, ptr, col, val), rhs, prm); \
    break;

  switch (cb) {
  case 1:
    solve_cpr(std::tie(rows, ptr, col, val), rhs, prm);
    break;

    BOOST_PP_SEQ_FOR_EACH(CALL_BLOCK_SOLVER, ~, ARCANE_ALINA_BLOCK_SIZES)

  default:
    precondition(false, "Unsupported block size");
    break;
  }

  std::cout << prof << std::endl;
}
