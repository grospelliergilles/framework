// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedRelaxationRuntime.h                              (C) 2026-2026 */
/*                                                                           */
/* Distributed memory sparse approximate inverse relaxation scheme.          */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MPI_DISTRIBUTEDRELAXATIONRUNTIME_H
#define ARCANE_ALINA_MPI_DISTRIBUTEDRELAXATIONRUNTIME_H
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

#include <arcane/alina/BackendInterface.h>
#include <arcane/alina/RelaxationRuntime.h>
#include <arcane/alina/DistributedRelaxation.h>
#include <arcane/alina/DistributedMatrix.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::mpi::relaxation
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Distributed memory sparse approximate inverse relaxation scheme.
 */
template <class Backend>
struct DistributedRelaxationRuntime
{
  typedef Backend backend_type;
  typedef typename Backend::params backend_params;
  typedef Alina::PropertyTree params;

  runtime::relaxation::eRelaxationType r;
  void* handle = nullptr;

  DistributedRelaxationRuntime(const DistributedMatrix<Backend>& A,
                               params prm, const backend_params& bprm = backend_params())
  : r(prm.get("type", runtime::relaxation::eRelaxationType::spai0))
  {
    if (!prm.erase("type"))
      ARCANE_ALINA_PARAM_MISSING("type");

    switch (r) {

#define ARCANE_ALINA_RELAX_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    handle = static_cast<void*>(new ::Arcane::Alina::mpi::relaxation::Distributed##type<Backend>(A, prm, bprm)); \
    break

#define ARCANE_ALINA_RELAX_LOCAL_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    handle = call_constructor<::Arcane::Alina::relaxation::type>(A, prm, bprm); \
    break;

#define ARCANE_ALINA_RELAX_LOCAL_LOCAL(type) \
  case runtime::relaxation::eRelaxationType::type: \
    handle = call_constructor<::Arcane::Alina::relaxation::type>(*A.local(), prm, bprm); \
    break;

      ARCANE_ALINA_RELAX_DISTR(SPAI0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ChebyshevRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(DampedJacobiRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(ILU0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(ILUKRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(ILUPRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(ILUTRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(SPAI1Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(GaussSeidelRelaxation);

#undef ARCANE_ALINA_RELAX_LOCAL_LOCAL
#undef ARCANE_ALINA_RELAX_LOCAL_DISTR
#undef ARCANE_ALINA_RELAX_DISTR

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  ~DistributedRelaxationRuntime()
  {
    switch (r) {
#define ARCANE_ALINA_RELAX_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    delete static_cast<::Arcane::Alina::mpi::relaxation::Distributed##type<Backend>*>(handle); \
    break

#define ARCANE_ALINA_RELAX_LOCAL(type) \
  case runtime::relaxation::eRelaxationType::type: \
    delete static_cast<::Arcane::Alina::relaxation::type<Backend>*>(handle); \
    break;

      ARCANE_ALINA_RELAX_DISTR(SPAI0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL(DampedJacobiRelaxation);
      ARCANE_ALINA_RELAX_LOCAL(ILU0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL(ILUKRelaxation);
      ARCANE_ALINA_RELAX_LOCAL(ILUPRelaxation);
      ARCANE_ALINA_RELAX_LOCAL(ILUTRelaxation);
      ARCANE_ALINA_RELAX_LOCAL(SPAI1Relaxation);
      ARCANE_ALINA_RELAX_LOCAL(ChebyshevRelaxation);
      ARCANE_ALINA_RELAX_LOCAL(GaussSeidelRelaxation);

#undef ARCANE_ALINA_RELAX_LOCAL
#undef ARCANE_ALINA_RELAX_DISTR

    default:
      break;
    }
  }

  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_pre(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    switch (r) {

#define ARCANE_ALINA_RELAX_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    static_cast<const ::Arcane::Alina::mpi::relaxation::Distributed##type<Backend>*>(handle)->apply_pre(A, rhs, x, tmp); \
    break

#define ARCANE_ALINA_RELAX_LOCAL_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    call_apply_pre<::Arcane::Alina::relaxation::type>(A, rhs, x, tmp); \
    break;

#define ARCANE_ALINA_RELAX_LOCAL_LOCAL(type) \
  case runtime::relaxation::eRelaxationType::type: \
    call_apply_pre<::Arcane::Alina::relaxation::type>(*A.local_backend(), rhs, x, tmp); \
    break;

      ARCANE_ALINA_RELAX_DISTR(SPAI0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(DampedJacobiRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILU0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUKRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUPRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUTRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(SPAI1Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ChebyshevRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(GaussSeidelRelaxation);

#undef ARCANE_ALINA_RELAX_LOCAL_LOCAL
#undef ARCANE_ALINA_RELAX_LOCAL_DISTR
#undef ARCANE_ALINA_RELAX_DISTR

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_post(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    switch (r) {

#define ARCANE_ALINA_RELAX_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    static_cast<const ::Arcane::Alina::mpi::relaxation::Distributed##type<Backend>*>(handle)->apply_post(A, rhs, x, tmp); \
    break

#define ARCANE_ALINA_RELAX_LOCAL_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    call_apply_post<::Arcane::Alina::relaxation::type>(A, rhs, x, tmp); \
    break;

#define ARCANE_ALINA_RELAX_LOCAL_LOCAL(type) \
  case runtime::relaxation::eRelaxationType::type: \
    call_apply_post<::Arcane::Alina::relaxation::type>(*A.local_backend(), rhs, x, tmp); \
    break;

      ARCANE_ALINA_RELAX_DISTR(SPAI0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(DampedJacobiRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILU0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUKRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUPRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUTRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(SPAI1Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ChebyshevRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(GaussSeidelRelaxation);

#undef ARCANE_ALINA_RELAX_LOCAL_LOCAL
#undef ARCANE_ALINA_RELAX_LOCAL_DISTR
#undef ARCANE_ALINA_RELAX_DISTR

    default:
      throw std::invalid_argument("Unsupported relaxation type");
    }
  }

  template <class Matrix, class VectorRHS, class VectorX>
  void apply(const Matrix& A, const VectorRHS& rhs, VectorX& x) const
  {
    switch (r) {

#define ARCANE_ALINA_RELAX_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    static_cast<const ::Arcane::Alina::mpi::relaxation::Distributed##type<Backend>*>(handle)->apply(A, rhs, x); \
    break

#define ARCANE_ALINA_RELAX_LOCAL_DISTR(type) \
  case runtime::relaxation::eRelaxationType::type: \
    call_apply<::Arcane::Alina::relaxation::type>(A, rhs, x); \
    break;

#define ARCANE_ALINA_RELAX_LOCAL_LOCAL(type) \
  case runtime::relaxation::eRelaxationType::type: \
    call_apply<::Arcane::Alina::relaxation::type>(*A.local_backend(), rhs, x); \
    break;

      ARCANE_ALINA_RELAX_DISTR(SPAI0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(DampedJacobiRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_LOCAL(GaussSeidelRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILU0Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUKRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUPRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ILUTRelaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(SPAI1Relaxation);
      ARCANE_ALINA_RELAX_LOCAL_DISTR(ChebyshevRelaxation);

#undef ARCANE_ALINA_RELAX_LOCAL_LOCAL
#undef ARCANE_ALINA_RELAX_LOCAL_DISTR
#undef ARCANE_ALINA_RELAX_DISTR

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

} // namespace Arcane::Alina::runtime::mpi::relaxation

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
