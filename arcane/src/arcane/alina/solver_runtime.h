#ifndef ARCANE_ALINA_SOLVER_RUNTIME_HPP
#define ARCANE_ALINA_SOLVER_RUNTIME_HPP

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

/**
 * \file   alina/solver/runtime.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Runtime-configurable wrappers around amgcl iterative solvers.
 */

#include <iostream>
#include <stdexcept>
#include <type_traits>

#include <arcane/alina/util.h>
#include <arcane/alina/ConjugateGradientSolver.h>
#include <arcane/alina/BiCGStabSolver.h>
#include <arcane/alina/BiCGStabLSolver.h>
#include <arcane/alina/GMRESSolver.h>
#include <arcane/alina/LooseGMRESSolver.h>
#include <arcane/alina/FlexibleGMRESSolver.h>
#include <arcane/alina/IDRSSolver.h>
#include <arcane/alina/RichardsonSolver.h>
#include <arcane/alina/PreconditionerOnlySolver.h>
#include <arcane/alina/solver_detail_default_inner_product.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::solver
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

enum type
{
  cg, ///< Conjugate gradients method
  ConjugateGradientSolver = cg, ///< Conjugate gradients method
  bicgstab, ///< BiConjugate Gradient Stabilized
  BiCGStabSolver = bicgstab, ///< BiConjugate Gradient Stabilized
  bicgstabl, ///< BiCGStab(ell)
  BiCGStabLSolver = bicgstabl, ///< BiCGStab(ell)
  gmres, ///< GMRES
  GMRESSolver = gmres, ///< GMRES
  lgmres, ///< LGMRES
  LooseGMRESSolver = lgmres, ///< LGMRES
  fgmres, ///< FGMRES
  FlexibleGMRESSolver = fgmres,
  idrs, ///< IDR(s)
  IDRSSolver = idrs, ///< IDR(s)
  richardson, ///< Richardson iteration
  RichardsonSolver = richardson, ///< Richardson iteration
  preonly, ///< Only apply preconditioner once
  PreconditionerOnlySolver = preonly
};

inline std::ostream& operator<<(std::ostream& os, type s)
{
  switch (s) {
  case cg:
    return os << "cg";
  case bicgstab:
    return os << "bicgstab";
  case bicgstabl:
    return os << "bicgstabl";
  case gmres:
    return os << "gmres";
  case lgmres:
    return os << "lgmres";
  case fgmres:
    return os << "fgmres";
  case idrs:
    return os << "idrs";
  case richardson:
    return os << "richardson";
  case preonly:
    return os << "preonly";
  default:
    return os << "???";
  }
}

inline std::istream& operator>>(std::istream& in, type& s)
{
  std::string val;
  in >> val;

  if (val == "cg")
    s = cg;
  else if (val == "bicgstab")
    s = bicgstab;
  else if (val == "bicgstabl")
    s = bicgstabl;
  else if (val == "gmres")
    s = gmres;
  else if (val == "lgmres")
    s = lgmres;
  else if (val == "fgmres")
    s = fgmres;
  else if (val == "idrs")
    s = idrs;
  else if (val == "richardson")
    s = richardson;
  else if (val == "preonly")
    s = preonly;
  else
    throw std::invalid_argument("Invalid solver value. Valid choices are: "
                                "cg, bicgstab, bicgstabl, gmres, lgmres, fgmres, idrs, richardson, preonly.");

  return in;
}

#define ARCANE_ALINA_ALL_RUNTIME_SOLVER() \
  ARCANE_ALINA_RUNTIME_SOLVER(ConjugateGradientSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(BiCGStabSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(BiCGStabLSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(GMRESSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(LooseGMRESSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(FlexibleGMRESSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(IDRSSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(RichardsonSolver); \
  ARCANE_ALINA_RUNTIME_SOLVER(PreconditionerOnlySolver)

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend,
          class InnerProduct = Alina::solver::detail::default_inner_product>
struct wrapper
{
  typedef PropertyTree params;
  typedef typename Backend::params backend_params;
  typedef typename Backend::value_type value_type;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  typedef Backend backend_type;

  type s;
  void* handle = nullptr;

  explicit wrapper(size_t n, params prm = params(),
                   const backend_params& bprm = backend_params(),
                   const InnerProduct& inner_product = InnerProduct())
  : s(prm.get("type", runtime::solver::bicgstab))
  {
    if (!prm.erase("type"))
      ARCANE_ALINA_PARAM_MISSING("type");

    switch (s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case type: \
    handle = static_cast<void*>(new ::Arcane::Alina::solver::type<Backend, InnerProduct>(n, prm, bprm, inner_product)); \
    break

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER

    default:
      throw std::invalid_argument("Unsupported solver type");
    }
  }

  ~wrapper()
  {
    switch (s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case type: \
    delete static_cast<Alina::solver::type<Backend, InnerProduct>*>(handle); \
    break

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER
    }
  }

  template <class Matrix, class Precond, class Vec1, class Vec2>
  SolverResult operator()(const Matrix& A, const Precond& P, const Vec1& rhs, Vec2&& x) const
  {
    switch (s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case type: \
    return static_cast<Alina::solver::type<Backend, InnerProduct>*>(handle)->operator()(A, P, rhs, x)

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER

    default:
      throw std::invalid_argument("Unsupported solver type");
    }
  }

  template <class Precond, class Vec1, class Vec2>
  SolverResult operator()(const Precond& P, const Vec1& rhs, Vec2&& x) const
  {
    return (*this)(P.system_matrix(), P, rhs, x);
  }

  friend std::ostream& operator<<(std::ostream& os, const wrapper& w)
  {
    switch (w.s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case type: \
    return os << *static_cast<Alina::solver::type<Backend, InnerProduct>*>(w.handle)

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER

    default:
      throw std::invalid_argument("Unsupported solver type");
    }
  }

  size_t bytes() const
  {
    switch (s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case type: \
    return backend::bytes(*static_cast<Alina::solver::type<Backend, InnerProduct>*>(handle))

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER

    default:
      throw std::invalid_argument("Unsupported solver type");
    }
  }
};

#undef ARCANE_ALINA_ALL_RUNTIME_SOLVER

} // namespace Arcane::Alina::runtime::solver

#endif
