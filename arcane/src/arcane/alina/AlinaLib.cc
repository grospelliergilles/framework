#include <iostream>

#include <type_traits>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <arcane/alina/RelaxationRuntime.h>
#include <arcane/alina/coarsening_runtime.h>
#include <arcane/alina/solver_runtime.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/AMG.h>
#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/Adapters.h>

#include "AlinaLib.h"

#ifdef ARCANE_ALINA_PROFILING
#include <arcane/alina/profiler.h>
namespace amgcl {
profiler<> prof;
}
#endif

using namespace Arcane;

//---------------------------------------------------------------------------

typedef Alina::backend::BuiltinBackend<double> Backend;
typedef Alina::AMG<Backend, Alina::runtime::coarsening::wrapper, Alina::runtime::relaxation::wrapper> AMG;
typedef Alina::runtime::solver::wrapper<Backend>  ISolver;
typedef Alina::make_solver<AMG, ISolver> Solver;
typedef Alina::PropertyTree Params;

conv_info _toConvInfo(const Alina::SolverResult& r)
{
  conv_info x;
  x.iterations = r.nbIteration();
  x.residual = r.residual();
  return x;
}

//---------------------------------------------------------------------------

amgclHandle STDCALL ARCANE_ALINA_params_create()
{
  return static_cast<amgclHandle>(new Params());
}

//---------------------------------------------------------------------------
void STDCALL ARCANE_ALINA_params_seti(amgclHandle prm, const char* name, int value)
{
  static_cast<Params*>(prm)->put(name, value);
}

//---------------------------------------------------------------------------
void STDCALL ARCANE_ALINA_params_setf(amgclHandle prm, const char* name, float value)
{
  static_cast<Params*>(prm)->put(name, value);
}

//---------------------------------------------------------------------------
void STDCALL ARCANE_ALINA_params_sets(amgclHandle prm, const char* name, const char* value)
{
  static_cast<Params*>(prm)->put(name, value);
}

//---------------------------------------------------------------------------
void STDCALL ARCANE_ALINA_params_read_json(amgclHandle prm, const char* fname)
{
  Params& p = *static_cast<Params*>(prm);
  p.read_json(fname);
}

//---------------------------------------------------------------------------
void STDCALL ARCANE_ALINA_params_destroy(amgclHandle prm)
{
  delete static_cast<Params*>(prm);
}

//---------------------------------------------------------------------------
amgclHandle STDCALL
ARCANE_ALINA_precond_create(int n,
                     const int* ptr,
                     const int* col,
                     const double* val,
                     amgclHandle prm)
{
  auto A = std::make_tuple(n,
                           boost::make_iterator_range(ptr, ptr + n + 1),
                           boost::make_iterator_range(col, col + ptr[n]),
                           boost::make_iterator_range(val, val + ptr[n]));

  if (prm)
    return static_cast<amgclHandle>(new AMG(A, *static_cast<Params*>(prm)));
  else
    return static_cast<amgclHandle>(new AMG(A));
}

//---------------------------------------------------------------------------
void STDCALL
ARCANE_ALINA_precond_apply(amgclHandle handle, const double* rhs, double* x)
{
  AMG* amg = static_cast<AMG*>(handle);

  size_t n = Alina::backend::rows(amg->system_matrix());

  boost::iterator_range<double*> x_range =
  boost::make_iterator_range(x, x + n);

  amg->apply(boost::make_iterator_range(rhs, rhs + n), x_range);
}

//---------------------------------------------------------------------------
void STDCALL
ARCANE_ALINA_precond_report(amgclHandle handle)
{
  std::cout << *static_cast<AMG*>(handle) << std::endl;
}

//---------------------------------------------------------------------------
void STDCALL
ARCANE_ALINA_precond_destroy(amgclHandle handle)
{
  delete static_cast<AMG*>(handle);
}

//---------------------------------------------------------------------------
amgclHandle STDCALL
ARCANE_ALINA_solver_create(int n, const int* ptr,
                           const int* col,
                           const double* val,
                           amgclHandle prm)
{
  auto A = std::make_tuple(n,
                           boost::make_iterator_range(ptr, ptr + n + 1),
                           boost::make_iterator_range(col, col + ptr[n]),
                           boost::make_iterator_range(val, val + ptr[n]));

  if (prm)
    return static_cast<amgclHandle>(new Solver(A, *static_cast<Params*>(prm)));
  else
    return static_cast<amgclHandle>(new Solver(A));
}

//---------------------------------------------------------------------------
void STDCALL
ARCANE_ALINA_solver_report(amgclHandle handle)
{
  std::cout << static_cast<Solver*>(handle)->precond() << std::endl;
}

//---------------------------------------------------------------------------
void STDCALL
ARCANE_ALINA_solver_destroy(amgclHandle handle)
{
  delete static_cast<Solver*>(handle);
}

//---------------------------------------------------------------------------
conv_info STDCALL
ARCANE_ALINA_solver_solve(amgclHandle handle,
                          const double* rhs,
                          double* x)
{
  Solver* slv = static_cast<Solver*>(handle);

  size_t n = slv->size();

  boost::iterator_range<double*> x_range = boost::make_iterator_range(x, x + n);

  Alina::SolverResult r = (*slv)(boost::make_iterator_range(rhs, rhs + n), x_range);

  return _toConvInfo(r);
}

//---------------------------------------------------------------------------
conv_info STDCALL
ARCANE_ALINA_solver_solve_mtx(amgclHandle handle,
                              int const* A_ptr,
                              int const* A_col,
                              double const* A_val,
                              const double* rhs,
                              double* x)
{
  Solver* slv = static_cast<Solver*>(handle);

  size_t n = slv->size();

  boost::iterator_range<double*> x_range = boost::make_iterator_range(x, x + n);

  Alina::SolverResult r = (*slv)(
  std::make_tuple(n,
                  boost::make_iterator_range(A_ptr, A_ptr + n + 1),
                  boost::make_iterator_range(A_col, A_col + A_ptr[n]),
                  boost::make_iterator_range(A_val, A_val + A_ptr[n])),
  boost::make_iterator_range(rhs, rhs + n), x_range);

  return _toConvInfo(r);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
