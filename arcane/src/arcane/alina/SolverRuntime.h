// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* SolverRuntime.h                                             (C) 2026-2026 */
/*                                                                           */
/* Runtime-configurable solvers.                                             */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_SOLVERRUNTIME_H
#define ARCANE_ALINA_SOLVERRUNTIME_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
 * This file is based on the work on AMGCL library (version march 2026)
 * which can be found at https://github.com/ddemidov/amgcl.
 *
 * Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>
 * SPDX-License-Identifier: MIT
 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

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

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

enum class eSolverType
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

inline std::ostream& operator<<(std::ostream& os, eSolverType s)
{
  switch (s) {
  case eSolverType::cg:
    return os << "cg";
  case eSolverType::bicgstab:
    return os << "bicgstab";
  case eSolverType::bicgstabl:
    return os << "bicgstabl";
  case eSolverType::gmres:
    return os << "gmres";
  case eSolverType::lgmres:
    return os << "lgmres";
  case eSolverType::fgmres:
    return os << "fgmres";
  case eSolverType::idrs:
    return os << "idrs";
  case eSolverType::richardson:
    return os << "richardson";
  case eSolverType::preonly:
    return os << "preonly";
  default:
    return os << "???";
  }
}

inline std::istream& operator>>(std::istream& in, eSolverType& s)
{
  std::string val;
  in >> val;

  if (val == "cg")
    s = eSolverType::cg;
  else if (val == "bicgstab")
    s = eSolverType::bicgstab;
  else if (val == "bicgstabl")
    s = eSolverType::bicgstabl;
  else if (val == "gmres")
    s = eSolverType::gmres;
  else if (val == "lgmres")
    s = eSolverType::lgmres;
  else if (val == "fgmres")
    s = eSolverType::fgmres;
  else if (val == "idrs")
    s = eSolverType::idrs;
  else if (val == "richardson")
    s = eSolverType::richardson;
  else if (val == "preonly")
    s = eSolverType::preonly;
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
/*!
 * \brief Runtime-configurable wrappers around iterative solvers.
 */
template <class Backend, class InnerProduct = Alina::solver::detail::default_inner_product>
struct SolverRuntime
{
  typedef PropertyTree params;
  typedef typename Backend::params backend_params;
  typedef typename Backend::value_type value_type;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  typedef Backend backend_type;

  eSolverType s;
  void* handle = nullptr;

  explicit SolverRuntime(size_t n, params prm = params(),
                         const backend_params& bprm = backend_params(),
                         const InnerProduct& inner_product = InnerProduct())
  : s(prm.get("type", eSolverType::bicgstab))
  {
    if (!prm.erase("type"))
      ARCANE_ALINA_PARAM_MISSING("type");

    switch (s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case eSolverType::type: \
    handle = static_cast<void*>(new ::Arcane::Alina::solver::type<Backend, InnerProduct>(n, prm, bprm, inner_product)); \
    break

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER

    default:
      throw std::invalid_argument("Unsupported solver type");
    }
  }

  ~SolverRuntime()
  {
    switch (s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case eSolverType::type: \
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
  case eSolverType::type: \
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

  friend std::ostream& operator<<(std::ostream& os, const SolverRuntime& w)
  {
    switch (w.s) {

#define ARCANE_ALINA_RUNTIME_SOLVER(type) \
  case eSolverType::type: \
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
  case eSolverType::type: \
    return backend::bytes(*static_cast<Alina::solver::type<Backend, InnerProduct>*>(handle))

      ARCANE_ALINA_ALL_RUNTIME_SOLVER();

#undef ARCANE_ALINA_RUNTIME_SOLVER

    default:
      throw std::invalid_argument("Unsupported solver type");
    }
  }
};

#undef ARCANE_ALINA_ALL_RUNTIME_SOLVER

} // namespace Arcane::Alina

#endif
