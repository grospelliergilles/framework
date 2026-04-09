#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <numeric>

#include <boost/program_options.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/scope_exit.hpp>

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/PreconditionerRuntime.h>
#include <arcane/alina/Adapters.h>
#include <arcane/alina/mpi/mp_make_solver.h>
#include <arcane/alina/mpi/mp_solver_runtime.h>
#include <arcane/alina/mpi/mp_block_preconditioner.h>
#include <arcane/alina/profiler.h>
#include <arcane/alina/util.h>

// Pour test compilation uniquement
#include <arcane/alina/mpi/mp_solver.h>
#include <arcane/alina/mpi/mp_solver_runtime.h>

// Pour test compilation uniquement
#include <arcane/alina/mpi/mp_relaxation_runtime.h>

#include "domain_partition.h"

using namespace Arcane;

namespace Arcane::Alina { profiler<> prof; }
using Alina::prof;
using Alina::precondition;

//---------------------------------------------------------------------------
struct renumbering {
    const domain_partition<2> &part;
    const std::vector<ptrdiff_t> &dom;

    renumbering(
            const domain_partition<2> &p,
            const std::vector<ptrdiff_t> &d
            ) : part(p), dom(d)
    {}

    ptrdiff_t operator()(ptrdiff_t i, ptrdiff_t j) const {
        boost::array<ptrdiff_t, 2> p = {{i, j}};
        std::pair<int,ptrdiff_t> v = part.index(p);
        return dom[v.first] + v.second;
    }
};

//---------------------------------------------------------------------------
template <template <class> class Precond, class Matrix>
std::tuple<size_t, double> solve(
        const Alina::mpi::communicator &comm,
        const Alina::PropertyTree &prm,
        const Matrix &A
        )
{
    typedef Alina::backend::BuiltinBackend<double> Backend;

    typedef Alina::mpi::make_solver<
        Alina::mpi::block_preconditioner< Precond<Backend> >,
        Alina::runtime::mpi::solver::wrapper<Backend>
        > Solver;

    const size_t n = Alina::backend::rows(A);

    std::vector<double> rhs(n, 1), x(n, 0);

    prof.tic("setup");
    Solver solve(comm, A, prm);
    prof.toc("setup");

    {
        auto t2 = prof.scoped_tic("solve");
        return solve(rhs, x);
    }
}

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    namespace po = boost::program_options;

    using Alina::prof;
    using std::vector;
    using std::string;

    po::options_description desc("Options");

    desc.add_options()
        ("help,h", "Show this help.")
        ("prm-file,P",
         po::value<string>(),
         "Parameter file in json format. "
        )
        (
         "prm,p",
         po::value< vector<string> >()->multitoken(),
         "Parameters specified as name=value pairs. "
         "May be provided multiple times. Examples:\n"
         "  -p solver.tol=1e-3\n"
         "  -p precond.coarse_enough=300"
        )
        (
         "size,n",
         po::value<int>()->default_value(1024),
         "The size of the Poisson problem to solve. "
         "Specified as number of grid nodes along each dimension of a unit square. "
         "The resulting system will have n*n unknowns. "
        )
        (
         "single-level,1",
         po::bool_switch()->default_value(false),
         "When specified, the AMG hierarchy is not constructed. "
         "Instead, the problem is solved using a single-level smoother as preconditioner. "
        )
        (
         "initial,x",
         po::value<double>()->default_value(0),
         "Value to use as initial approximation. "
        )
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    Alina::PropertyTree prm;
    if (vm.count("prm-file")) {
      prm.read_json(vm["prm-file"].as<string>());
    }

    if (vm.count("prm")) {
        for(const string &v : vm["prm"].as<vector<string> >()) {
            Alina::put(prm, v);
        }
    }

    MPI_Init(&argc, &argv);
    BOOST_SCOPE_EXIT(void) {
        MPI_Finalize();
    } BOOST_SCOPE_EXIT_END

    Alina::mpi::communicator world(MPI_COMM_WORLD);

    if (world.rank == 0)
        std::cout << "World size: " << world.size << std::endl;

    const ptrdiff_t n   = vm["size"].as<int>();
    const double    h2i = (n - 1) * (n - 1);

    boost::array<ptrdiff_t, 2> lo = { {0, 0} };
    boost::array<ptrdiff_t, 2> hi = { {n - 1, n - 1} };

    prof.tic("partition");
    domain_partition<2> part(lo, hi, world.size);
    ptrdiff_t chunk = part.size( world.rank );

    std::vector<ptrdiff_t> domain(world.size + 1);
    MPI_Allgather(
            &chunk, 1, Alina::mpi::datatype<ptrdiff_t>(),
            &domain[1], 1, Alina::mpi::datatype<ptrdiff_t>(), world);
    std::partial_sum(domain.begin(), domain.end(), domain.begin());

    lo = part.domain(world.rank).min_corner();
    hi = part.domain(world.rank).max_corner();
    prof.toc("partition");

    renumbering renum(part, domain);

    prof.tic("assemble");
    std::vector<ptrdiff_t> ptr;
    std::vector<ptrdiff_t> col;
    std::vector<double>    val;
    std::vector<double>    rhs;

    ptr.reserve(chunk + 1);
    col.reserve(chunk * 5);
    val.reserve(chunk * 5);

    ptr.push_back(0);

    for(ptrdiff_t j = lo[1]; j <= hi[1]; ++j) {
        for(ptrdiff_t i = lo[0]; i <= hi[0]; ++i) {
            if (j > 0)  {
                col.push_back(renum(i,j-1));
                val.push_back(-h2i);
            }

            if (i > 0) {
                col.push_back(renum(i-1,j));
                val.push_back(-h2i);
            }

            col.push_back(renum(i,j));
            val.push_back(4 * h2i);

            if (i + 1 < n) {
                col.push_back(renum(i+1,j));
                val.push_back(-h2i);
            }

            if (j + 1 < n) {
                col.push_back(renum(i,j+1));
                val.push_back(-h2i);
            }

            ptr.push_back( col.size() );
        }
    }
    prof.toc("assemble");

    size_t iters;
    double error;

    bool single_level = vm["single-level"].as<bool>();

    if (single_level)
        prm.put("precond.class", "relaxation");

    std::tie(iters, error) = solve<Alina::runtime::PreconditionerRuntime>(
            world, prm, std::tie(chunk, ptr, col, val));

    if (world.rank == 0) {
        std::cout
            << "Iterations: " << iters << std::endl
            << "Error:      " << error << std::endl
            << std::endl
            << prof << std::endl;
    }
}
