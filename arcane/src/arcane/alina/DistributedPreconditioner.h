// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedPreconditioner.h                                 (C) 2026-2026 */
/*                                                                           */
/* Runtime wrapper around mpi preconditioners.                               */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MPI_DISTRIBUTEDPRECONDITIONER_H
#define ARCANE_ALINA_MPI_DISTRIBUTEDPRECONDITIONER_H
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

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/DistributedAMG.h>
#include <arcane/alina/DistributedInnerProduct.h>
#include <arcane/alina/DistributedCoarseningRuntime.h>
#include <arcane/alina/DistributedRelaxationRuntime.h>
#include <arcane/alina/DistributedDirectSolverRuntime.h>
#include <arcane/alina/MatrixPartitionerRuntime.h>
#include <arcane/alina/DistributedRelaxation.h>
#include <arcane/alina/DistributedMatrix.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::mpi
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Preconditioner kinds.
namespace precond_class
{
  enum type
  {
    amg, ///< AMG
    relaxation ///< Single-level relaxation
  };

  inline std::ostream& operator<<(std::ostream& os, type p)
  {
    switch (p) {
    case amg:
      return os << "amg";
    case relaxation:
      return os << "relaxation";
    default:
      return os << "???";
    }
  }

  inline std::istream& operator>>(std::istream& in, type& p)
  {
    std::string val;
    in >> val;

    if (val == "amg")
      p = amg;
    else if (val == "relaxation")
      p = relaxation;
    else
      throw std::invalid_argument("Invalid preconditioner class. "
                                  "Valid choices are: amg, relaxation");

    return in;
  }
} // namespace precond_class

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Distributed Preconditioner.
 */
template <class Backend>
class DistributedPreconditioner
{
 public:

  typedef Backend backend_type;
  typedef typename backend_type::params backend_params;
  typedef PropertyTree params;
  typedef typename backend_type::value_type value_type;
  typedef DistributedMatrix<backend_type> matrix;

  using AMGPrecondType = Alina::mpi::DistributedAMG<Backend,
                                                    Alina::runtime::mpi::coarsening::DistributedCoarseningRuntime<Backend>,
                                                    Alina::runtime::mpi::relaxation::DistributedRelaxationRuntime<Backend>,
                                                    DistributedDirectSolverRuntime<value_type>,
                                                    MatrixPartitionerRuntime<Backend>>;

  template <class Matrix>
  DistributedPreconditioner(mpi_communicator comm,
                            const Matrix& Astrip,
                            params prm = params(),
                            const backend_params& bprm = backend_params())
  : _class(prm.get("class", precond_class::amg))
  , handle(0)
  {
    init(std::make_shared<matrix>(comm, Astrip, backend::rows(Astrip)), prm, bprm);
  }

  DistributedPreconditioner(mpi_communicator,
                            std::shared_ptr<matrix> A,
                            params prm = params(),
                            const backend_params& bprm = backend_params())
  : _class(prm.get("class", precond_class::amg))
  , handle(0)
  {
    init(A, prm, bprm);
  }

  ~DistributedPreconditioner()
  {
    switch (_class) {
    case precond_class::amg: {
      delete static_cast<AMGPrecondType*>(handle);
    } break;
    case precond_class::relaxation: {
      typedef Alina::mpi::relaxation::as_preconditioner<
      Alina::runtime::mpi::relaxation::DistributedRelaxationRuntime<Backend>>
      Precond;

      delete static_cast<Precond*>(handle);
    } break;
    default:
      break;
    }
  }

  template <class Matrix>
  void rebuild(const Matrix& A,
               const backend_params& bprm = backend_params())
  {
    switch (_class) {
    case precond_class::amg: {
      static_cast<AMGPrecondType*>(handle)->rebuild(A, bprm);
    } break;
    default:
      std::cerr << "rebuild is a noop unless the preconditioner is AMG" << std::endl;
      return;
    }
  }

  template <class Vec1, class Vec2>
  void apply(const Vec1& rhs, Vec2&& x) const
  {
    switch (_class) {
    case precond_class::amg: {
      static_cast<AMGPrecondType*>(handle)->apply(rhs, x);
    } break;
    case precond_class::relaxation: {
      typedef Alina::mpi::relaxation::as_preconditioner<
      Alina::runtime::mpi::relaxation::DistributedRelaxationRuntime<Backend>>
      Precond;

      static_cast<Precond*>(handle)->apply(rhs, x);
    } break;
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

  /// Returns the system matrix from the finest level.
  std::shared_ptr<matrix> system_matrix_ptr() const
  {
    switch (_class) {
    case precond_class::amg: {
      return static_cast<AMGPrecondType*>(handle)->system_matrix_ptr();
    }
    case precond_class::relaxation: {
      typedef Alina::mpi::relaxation::as_preconditioner<
      Alina::runtime::mpi::relaxation::DistributedRelaxationRuntime<Backend>>
      Precond;

      return static_cast<Precond*>(handle)->system_matrix_ptr();
    }
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

  const matrix& system_matrix() const
  {
    return *system_matrix_ptr();
  }

  friend std::ostream& operator<<(std::ostream& os, const DistributedPreconditioner& p)
  {
    switch (p._class) {
    case precond_class::amg: {
      return os << *static_cast<AMGPrecondType*>(p.handle);
    }
    case precond_class::relaxation: {
      typedef Alina::mpi::relaxation::as_preconditioner<
      Alina::runtime::mpi::relaxation::DistributedRelaxationRuntime<Backend>>
      Precond;

      return os << *static_cast<Precond*>(p.handle);
    }
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

 private:

  precond_class::type _class;
  void* handle;

  void init(std::shared_ptr<matrix> A, params& prm, const backend_params& bprm)
  {
    if (!prm.erase("class"))
      ARCANE_ALINA_PARAM_MISSING("class");

    switch (_class) {
    case precond_class::amg: {
      handle = static_cast<void*>(new AMGPrecondType(A->comm(), A, prm, bprm));
    } break;
    case precond_class::relaxation: {
      typedef Alina::mpi::relaxation::as_preconditioner<
      Alina::runtime::mpi::relaxation::DistributedRelaxationRuntime<Backend>>
      Precond;

      handle = static_cast<void*>(new Precond(A->comm(), A, prm, bprm));
    } break;
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::runtime::mpi

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::mpi
{

/*!
 * \brief Distributed block preconditioner.
 */
template <class Precond>
class block_preconditioner
{
 public:

  typedef typename Precond::params params;
  typedef typename Precond::backend_type backend_type;
  typedef typename backend_type::params backend_params;

  typedef typename backend_type::value_type value_type;
  typedef typename backend_type::matrix bmatrix;
  typedef DistributedMatrix<backend_type> matrix;

  template <class Matrix>
  block_preconditioner(mpi_communicator comm,
                       const Matrix& Astrip,
                       const params& prm = params(),
                       const backend_params& bprm = backend_params())
  {
    A = std::make_shared<matrix>(comm, Astrip, backend::rows(Astrip));
    P = std::make_shared<Precond>(A->local(), prm, bprm);
    A->set_local(P->system_matrix_ptr());
    A->move_to_backend(bprm);
  }

  block_preconditioner(mpi_communicator,
                       std::shared_ptr<matrix> A,
                       const params& prm = params(),
                       const backend_params& bprm = backend_params())
  : A(A)
  {
    P = std::make_shared<Precond>(A->local(), prm, bprm);
    A->set_local(P->system_matrix_ptr());
    A->move_to_backend(bprm);
  }

  std::shared_ptr<matrix> system_matrix_ptr() const
  {
    return A;
  }

  const matrix& system_matrix() const
  {
    return *A;
  }

  template <class Vec1, class Vec2>
  void apply(const Vec1& rhs, Vec2&& x) const
  {
    P->apply(rhs, x);
  }

 private:

  std::shared_ptr<matrix> A;
  std::shared_ptr<Precond> P;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::mpi

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
