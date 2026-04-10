#include <iostream>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/program_options.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/scope_exit.hpp>

#if defined(SOLVER_BACKEND_CUDA)
#  include <arcane/alina/CudaBackend.h>
#  include <arcane/alina/relaxation_cusparse_ilu0.h>
typedef Arcane::Alina::backend::cuda<double> Backend;
#else
#  ifndef SOLVER_BACKEND_BUILTIN
#    define SOLVER_BACKEND_BUILTIN
#  endif
#include <arcane/alina/BuiltinBackend.h>
typedef Arcane::Alina::backend::BuiltinBackend<double> Backend;
#endif

#include <arcane/alina/IO.h>
#include <arcane/alina/Adapters.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/CoarseningRuntime.h>
#include <arcane/alina/RelaxationRuntime.h>
#include <arcane/alina/mpi/mp_make_solver.h>
#include <arcane/alina/mpi/DistributedSchurPressureCorrection.h>
#include <arcane/alina/mpi/DistributedPreconditioner.h>
#include <arcane/alina/mpi/mp_subdomain_deflation.h>
#include <arcane/alina/mpi/DistributedSolverRuntime.h>
#include <arcane/alina/mpi/mp_direct_solver_runtime.h>
#include <arcane/alina/profiler.h>

namespace Arcane::Alina {
    profiler<> prof;
}

using namespace Arcane;

using Alina::prof;
using Alina::precondition;

//---------------------------------------------------------------------------
std::vector<ptrdiff_t> read_problem(
        const Alina::mpi::communicator &world,
        const std::string &A_file,
        const std::string &rhs_file,
        const std::string &part_file,
        std::vector<ptrdiff_t> &ptr,
        std::vector<ptrdiff_t> &col,
        std::vector<double>    &val,
        std::vector<double>    &rhs
        )
{
    // Read partition
    ptrdiff_t n, m;
    std::vector<ptrdiff_t> domain(world.size + 1, 0);
    std::vector<int> part;

    std::tie(n, m) = Alina::IO::mm_reader(part_file)(part);
    for(int p : part) {
        ++domain[p+1];
        precondition(p < world.size, "MPI world does not correspond to partition");
    }
    std::partial_sum(domain.begin(), domain.end(), domain.begin());

    ptrdiff_t chunk_beg = domain[world.rank];
    ptrdiff_t chunk_end = domain[world.rank + 1];
    ptrdiff_t chunk     = chunk_end - chunk_beg;

    // Reorder unknowns
    std::vector<ptrdiff_t> order(n);
    for(ptrdiff_t i = 0; i < n; ++i)
        order[i] = domain[part[i]]++;

    std::rotate(domain.begin(), domain.end()-1, domain.end());
    domain[0] = 0;

    // Read matrix chunk
    {
        using namespace Arcane::Alina::IO;

        std::ifstream A(A_file.c_str(), std::ios::binary);
        precondition(A, "Failed to open matrix file (" + A_file + ")");

        std::ifstream b(rhs_file.c_str(), std::ios::binary);
        precondition(b, "Failed to open rhs file (" + rhs_file + ")");

        ptrdiff_t rows;
        precondition(read(A, rows), "File I/O error");
        precondition(rows == n, "Matrix and partition have incompatible sizes");

        ptr.clear(); ptr.reserve(chunk + 1); ptr.push_back(0);

        std::vector<ptrdiff_t> gptr(n + 1);
        precondition(read(A, gptr), "File I/O error");

        size_t col_beg = sizeof(rows) + sizeof(gptr[0]) * (n + 1);
        size_t val_beg = col_beg + sizeof(col[0]) * gptr.back();
        size_t rhs_beg = 2 * sizeof(ptrdiff_t);

        // Count local nonzeros
        for(ptrdiff_t i = 0; i < n; ++i)
            if (part[i] == world.rank)
                ptr.push_back(gptr[i+1] - gptr[i]);

        std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());

        col.clear(); col.reserve(ptr.back());
        val.clear(); val.reserve(ptr.back());
        rhs.clear(); rhs.reserve(chunk);

        // Read local matrix and rhs stripes
        for(ptrdiff_t i = 0; i < n; ++i) {
            if (part[i] != world.rank) continue;

            ptrdiff_t c;
            A.seekg(col_beg + gptr[i] * sizeof(c));
            for(ptrdiff_t j = gptr[i], e = gptr[i+1]; j < e; ++j) {
                precondition(read(A, c), "File I/O error (1)");
                col.push_back(order[c]);
            }
        }

        for(ptrdiff_t i = 0; i < n; ++i) {
            if (part[i] != world.rank) continue;

            double v;
            A.seekg(val_beg + gptr[i] * sizeof(v));
            for(ptrdiff_t j = gptr[i], e = gptr[i+1]; j < e; ++j) {
                precondition(read(A, v), "File I/O error (2)");
                val.push_back(v);
            }
        }

        for(ptrdiff_t i = 0; i < n; ++i) {
            if (part[i] != world.rank) continue;

            double f;
            b.seekg(rhs_beg + i * sizeof(f));
            precondition(read(b, f), "File I/O error (3)");
            rhs.push_back(f);
        }
    }

    return domain;
}

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    BOOST_SCOPE_EXIT(void) {
        MPI_Finalize();
    } BOOST_SCOPE_EXIT_END

    Alina::mpi::communicator world(MPI_COMM_WORLD);

    if (world.rank == 0)
        std::cout << "World size: " << world.size << std::endl;

    // Read configuration from command line
    namespace po = boost::program_options;
    using std::string;
    po::options_description desc("Options");

    desc.add_options()
        ("help,h", "show help")
        (
         "matrix,A",
         po::value<string>()->required(),
         "The system matrix in binary format"
        )
        (
         "rhs,f",
         po::value<string>(),
         "The right-hand side in binary format"
        )
        (
         "part,s",
         po::value<string>()->required(),
         "Partitioning of the problem in MatrixMarket format"
        )
        (
         "pmask,m",
         po::value<string>(),
         "The pressure mask in binary format. Or, if the parameter has "
         "the form '%n:m', then each (n+i*m)-th variable is treated as pressure."
        )
        (
         "params,P",
         po::value<string>(),
         "parameter file in json format"
        )
        (
         "prm,p",
         po::value< std::vector<string> >()->multitoken(),
         "Parameters specified as name=value pairs. "
         "May be provided multiple times. Examples:\n"
         "  -p solver.tol=1e-3\n"
         "  -p precond.coarse_enough=300"
        )
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        if (world.rank == 0)
            std::cout << desc << std::endl;
        return 0;
    }

    po::notify(vm);

    Alina::PropertyTree prm;
    if (vm.count("params"))
      prm.read_json(vm["params"].as<string>());

    if (vm.count("prm")) {
        for(const string &v : vm["prm"].as<std::vector<string> >()) {
            Alina::put(prm, v);
        }
    }

    prof.tic("read problem");
    std::vector<ptrdiff_t> ptr;
    std::vector<ptrdiff_t> col;
    std::vector<double>    val;
    std::vector<double>    rhs;

    std::vector<ptrdiff_t> domain = read_problem(
            world,
            vm["matrix"].as<string>(), vm["rhs"].as<string>(), vm["part"].as<string>(),
            ptr, col, val, rhs
            );

    ptrdiff_t chunk = domain[world.rank + 1] - domain[world.rank];
    prof.toc("read problem");

    std::vector<char> pm;
    if(vm.count("pmask")) {
        std::string pmask = vm["pmask"].as<string>();
        prm.put("precond.pmask_size", chunk);

        switch (pmask[0]) {
            case '%':
            case '<':
            case '>':
                prm.put("precond.pmask_pattern", pmask);
                break;
            default:
                precondition(false, "Pressure mask may only be set with a pattern");
        }
    }

    std::function<double(ptrdiff_t,unsigned)> dv = Alina::mpi::constant_deflation(1);
    prm.put("precond.psolver.num_def_vec", 1);
    prm.put("precond.psolver.def_vec", &dv);

    Backend::params bprm;

#if defined(SOLVER_BACKEND_VEXCL)
    vex::Context ctx(vex::Filter::Env);
    std::cout << ctx << std::endl;
    bprm.q = ctx;
#elif defined(SOLVER_BACKEND_CUDA)
    cusparseCreate(&bprm.cusparse_handle);
#endif

    auto f = Backend::copy_vector(rhs, bprm);
    auto x = Backend::create_vector(chunk, bprm);

    Alina::backend::clear(*x);

    MPI_Barrier(world);

    prof.tic("setup");
    typedef
        Alina::mpi::DistributedPreconditionedSolver<
            Alina::mpi::DistributedSchurPressureCorrection<
                Alina::mpi::DistributedPreconditionedSolver<
                    Alina::mpi::block_preconditioner<
                        Alina::relaxation::as_preconditioner<Backend, Alina::runtime::relaxation::RuntimeRelaxation>
                        >,
                    Alina::runtime::mpi::solver::DistributedSolverRuntime<Backend>
                    >,
                Alina::mpi::subdomain_deflation<
                    Alina::AMG<Backend, Alina::runtime::coarsening::CoarseningRuntime, Alina::runtime::relaxation::RuntimeRelaxation>,
                    Alina::runtime::mpi::solver::DistributedSolverRuntime<Backend>,
                    Alina::runtime::mpi::direct::solver<double>
                    >
                >,
            Alina::runtime::mpi::solver::DistributedSolverRuntime<Backend>
            > Solver;

    Solver solve(world, std::tie(chunk, ptr, col, val), prm, bprm);
    double tm_setup = prof.toc("setup");

    prof.tic("solve");
    Alina::SolverResult r = solve(*f, *x);
    double tm_solve = prof.toc("solve");

    if (world.rank == 0) {
      std::cout << "Iters: " << r.nbIteration() << std::endl
                << "Error: " << r.residual() << std::endl
                << prof << std::endl;

#ifdef _OPENMP
        int nt = omp_get_max_threads();
#else
        int nt = 1;
#endif
        std::ostringstream log_name;
        log_name << "schur_" << domain.back() << "_" << nt << "_" << world.size << ".txt";
        std::ofstream log(log_name.str().c_str(), std::ios::app);
        log << domain.back() << "\t" << nt << "\t" << world.size
            << "\t" << tm_setup << "\t" << tm_solve
            << "\t" << r.nbIteration() << "\t" << std::endl;
    }
}
