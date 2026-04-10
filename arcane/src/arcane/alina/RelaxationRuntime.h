// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* relaxation_runtime.h                                        (C) 2026-2026 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_RELAXATION_RUNTIME_H
#define ARCANE_ALINA_RELAXATION_RUNTIME_H
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

#include "arcane/alina/AlinaGlobal.h"

#include <type_traits>

#include <arcane/alina/util.h>
#include <arcane/alina/Relaxation.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::relaxation
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Relaxation schemes.
enum type
{
  gauss_seidel, ///< Gauss-Seidel smoothing
  GaussSeidelRelaxation = gauss_seidel,
  ilu0, ///< Incomplete LU with zero fill-in
  ILU0Relaxation = ilu0,
  iluk, ///< Level-based incomplete LU
  ILUKRelaxation = iluk,
  ilup, ///< Level-based incomplete LU (fill-in is determined from A^p pattern)
  ILUPRelaxation = ilup,
  ilut, ///< Incomplete LU with thresholding
  ILUTRelaxation = ilut,
  damped_jacobi, ///< Damped Jacobi
  DampedJacobiRelaxation = damped_jacobi,
  spai0, ///< Sparse approximate inverse of 0th order
  SPAI0Relaxation = spai0,
  spai1, ///< Sparse approximate inverse of 1st order
  SPAI1Relaxation = spai1,
  chebyshev, ///< Chebyshev relaxation
  ChebyshevRelaxation = chebyshev
};

extern "C++" ARCANE_ALINA_EXPORT
std::ostream& operator<<(std::ostream& os, type r);

extern "C++" ARCANE_ALINA_EXPORT
std::istream& operator>>(std::istream& in, type& r);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct RuntimeRelaxation
{
  typedef Alina::PropertyTree params;
  typedef typename Backend::params backend_params;
  type r;
  void* handle;

  template <class Matrix>
  RuntimeRelaxation(const Matrix& A, params prm = params(),
          const backend_params& bprm = backend_params())
  : r(prm.get("type", runtime::relaxation::spai0))
  , handle(0)
  {
    if (!prm.erase("type"))
      ARCANE_ALINA_PARAM_MISSING("type");
    switch (r) {

#define ARCANE_ALINA_RUNTIME_RELAXATION(type) \
  case type: \
    handle = call_constructor<::Arcane::Alina::relaxation::type>(A, prm, bprm); \
    break

      ARCANE_ALINA_RUNTIME_RELAXATION(GaussSeidelRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILU0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUKRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUPRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUTRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(DampedJacobiRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI1Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ChebyshevRelaxation);

#undef ARCANE_ALINA_RUNTIME_RELAXATION

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  ~RuntimeRelaxation()
  {
    switch (r) {

#define ARCANE_ALINA_RUNTIME_RELAXATION(type) \
  case type: \
    delete static_cast<::Arcane::Alina::relaxation::type<Backend>*>(handle); \
    break

      ARCANE_ALINA_RUNTIME_RELAXATION(GaussSeidelRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILU0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUKRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUPRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUTRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(DampedJacobiRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI1Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ChebyshevRelaxation);

#undef ARCANE_ALINA_RUNTIME_RELAXATION
    }
  }

  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_pre(
  const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    std::cout << "PreconditionerPreRelaxationType=" << r << "\n";

    switch (r) {

#define ARCANE_ALINA_RUNTIME_RELAXATION(type) \
  case type: \
    call_apply_pre<::Arcane::Alina::relaxation::type>(A, rhs, x, tmp); \
    break

      ARCANE_ALINA_RUNTIME_RELAXATION(GaussSeidelRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILU0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUKRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUPRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUTRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(DampedJacobiRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI1Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ChebyshevRelaxation);

#undef ARCANE_ALINA_RUNTIME_RELAXATION

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_post(
  const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    std::cout << "PreconditionerPostRelaxationType=" << r << "\n";
    switch (r) {

#define ARCANE_ALINA_RUNTIME_RELAXATION(type) \
  case type: \
    call_apply_post<::Arcane::Alina::relaxation::type>(A, rhs, x, tmp); \
    break

      ARCANE_ALINA_RUNTIME_RELAXATION(GaussSeidelRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILU0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUKRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUPRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUTRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(DampedJacobiRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI1Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ChebyshevRelaxation);

#undef ARCANE_ALINA_RUNTIME_RELAXATION

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  template <class Matrix, class VectorRHS, class VectorX>
  void apply(const Matrix& A, const VectorRHS& rhs, VectorX& x) const
  {
    std::cout << "PreconditionerRelaxationType=" << r << "\n";

    switch (r) {

#define ARCANE_ALINA_RUNTIME_RELAXATION(type) \
  case type: \
    call_apply<Arcane::Alina::relaxation::type>(A, rhs, x); \
    break

      ARCANE_ALINA_RUNTIME_RELAXATION(GaussSeidelRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILU0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUKRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUPRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUTRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(DampedJacobiRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI1Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ChebyshevRelaxation);

#undef ARCANE_ALINA_RUNTIME_RELAXATION

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  size_t bytes() const
  {
    switch (r) {

#define ARCANE_ALINA_RUNTIME_RELAXATION(type) \
  case type: \
    return backend::bytes(*static_cast<::Arcane::Alina::relaxation::type<Backend>*>(handle))

      ARCANE_ALINA_RUNTIME_RELAXATION(GaussSeidelRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILU0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUKRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUPRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ILUTRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(DampedJacobiRelaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI0Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(SPAI1Relaxation);
      ARCANE_ALINA_RUNTIME_RELAXATION(ChebyshevRelaxation);

#undef ARCANE_ALINA_RUNTIME_RELAXATION

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  template <template <class> class Relaxation, class Matrix>
  typename std::enable_if<backend::relaxation_is_supported<Backend, Relaxation>::value, void*>::type
  call_constructor(const Matrix& A, const params& prm, const backend_params& bprm)
  {
    return static_cast<void*>(new Relaxation<Backend>(A, prm, bprm));
  }

  template <template <class> class Relaxation, class Matrix>
  typename std::enable_if<!backend::relaxation_is_supported<Backend, Relaxation>::value, void*>::type
  call_constructor(const Matrix&, const params&, const backend_params&)
  {
    throw std::logic_error("The relaxation is not supported by the backend");
  }

  template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  typename std::enable_if<backend::relaxation_is_supported<Backend, Relaxation>::value, void>::type
  call_apply_pre(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    static_cast<Relaxation<Backend>*>(handle)->apply_pre(A, rhs, x, tmp);
  }

  template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  typename std::enable_if<!backend::relaxation_is_supported<Backend, Relaxation>::value, void>::type
  call_apply_pre(const Matrix&, const VectorRHS&, VectorX&, VectorTMP&) const
  {
    throw std::logic_error("The relaxation is not supported by the backend");
  }

  template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  typename std::enable_if<backend::relaxation_is_supported<Backend, Relaxation>::value, void>::type
  call_apply_post(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    static_cast<Relaxation<Backend>*>(handle)->apply_post(A, rhs, x, tmp);
  }

  template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  typename std::enable_if<!backend::relaxation_is_supported<Backend, Relaxation>::value, void>::type
  call_apply_post(const Matrix&, const VectorRHS&, VectorX&, VectorTMP&) const
  {
    throw std::logic_error("The relaxation is not supported by the backend");
  }

  template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX>
  typename std::enable_if<backend::relaxation_is_supported<Backend, Relaxation>::value, void>::type
  call_apply(const Matrix& A, const VectorRHS& rhs, VectorX& x) const
  {
    static_cast<Relaxation<Backend>*>(handle)->apply(A, rhs, x);
  }

  template <template <class> class Relaxation, class Matrix, class VectorRHS, class VectorX>
  typename std::enable_if<!backend::relaxation_is_supported<Backend, Relaxation>::value, void>::type
  call_apply(const Matrix&, const VectorRHS&, VectorX&) const
  {
    throw std::logic_error("The relaxation is not supported by the backend");
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::runtime::relaxation

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
