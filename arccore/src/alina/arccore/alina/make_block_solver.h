#ifndef ARCANE_ALINA_MAKE_BLOCK_SOLVER_HPP
#define ARCANE_ALINA_MAKE_BLOCK_SOLVER_HPP

#include "arccore/alina/BackendInterface.h"
#include "arccore/alina/Adapters.h"
#include "arccore/alina/StaticMatrix.h"
#include "arccore/alina/PreconditionedSolver.h"
#include "arccore/alina/AlinaUtils.h"

namespace Arcane::Alina
{

/* Creates solver that operates in non-scalar domain but may take scalar inputs
 * for the system matrix and the rhs/solution vectors.
 */
template <class Precond, class IterativeSolver>
class make_block_solver
{
 public:

  typedef typename Precond::backend_type backend_type;
  typedef typename backend_type::value_type value_type;
  typedef typename backend_type::params backend_params;
  typedef typename backend_type::vector vector;
  typedef typename math::scalar_of<value_type>::type scalar_type;

  typedef typename PreconditionedSolver<Precond, IterativeSolver>::params params;

  template <class Matrix>
  make_block_solver(const Matrix& A,
                    const params& prm = params(),
                    const backend_params& bprm = backend_params())
  {
    S = std::make_shared<Solver>(adapter::block_matrix<value_type>(A), prm, bprm);
  }

  template <class Matrix, class Vec1, class Vec2>
  SolverResult operator()(const Matrix& A, const Vec1& rhs, Vec2&& x) const
  {
    auto F = backend::reinterpret_as_rhs<value_type>(rhs);
    auto X = backend::reinterpret_as_rhs<value_type>(x);

    return (*S)(A, F, X);
  }

  template <class Vec1, class Vec2>
  SolverResult operator()(const Vec1& rhs, Vec2&& x) const
  {
    auto F = backend::reinterpret_as_rhs<value_type>(rhs);
    auto X = backend::reinterpret_as_rhs<value_type>(x);

    return (*S)(F, X);
  }

  std::shared_ptr<typename Precond::matrix> system_matrix_ptr() const
  {
    return S->system_matrix_ptr();
  }

  typename Precond::matrix const& system_matrix() const
  {
    return S->system_matrix();
  }

  friend std::ostream& operator<<(std::ostream& os, const make_block_solver& p)
  {
    return os << *p.S << std::endl;
  }

  size_t bytes() const
  {
    return backend::bytes(*S);
  }

 private:

  typedef PreconditionedSolver<Precond, IterativeSolver> Solver;
  std::shared_ptr<Solver> S;
};

} // namespace Arcane::Alina

#endif
