// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedEigenSparseLUDirectSolver.h                      (C) 2026-2026 */
/*                                                                           */
/* Distributed wrapper for Eigen::SparseLU solver.                           */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_DISTRIBUTEDEIGENSPARSELUDIRECTSOLVER_H
#define ARCANE_ALINA_DISTRIBUTEDEIGENSPARSELUDIRECTSOLVER_H
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

#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/solver_eigen.h>
#include <arcane/alina/mpi/mp_util.h>
#include <arcane/alina/mpi/DistributedDirectSolverBase.h>

#include <Eigen/SparseLU>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::mpi::direct
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Distributed wrapper for Eigen::SparseLU solver.
 *
 * This is a wrapper around Eigen SparseLU solver that provides a
 * distributed direct solver interface but always works sequentially.
 */
template <typename value_type>
class DistributedEigenSparseLUDirectSolver
: public DistributedDirectSolverBase<value_type, DistributedEigenSparseLUDirectSolver<value_type>>
{
 public:

  using EigenMatrix = Eigen::SparseMatrix<value_type, Eigen::ColMajor, int>;
  using Solver = solver::EigenSolver<Eigen::SparseLU<EigenMatrix>>;
  typedef typename Solver::params params;
  typedef backend::CSRMatrix<value_type> build_matrix;

  /// Constructor.
  template <class Matrix>
  DistributedEigenSparseLUDirectSolver(mpi_communicator comm, const Matrix& A,
                                       const params& prm = params())
  : prm(prm)
  {
    static_cast<Base*>(this)->init(comm, A);
  }

  static size_t coarse_enough()
  {
    return Base::coarse_enough();
  }

  int comm_size(int /*n*/) const
  {
    return 1;
  }

  void init(mpi_communicator, const build_matrix& A)
  {
    S = std::make_shared<Solver>(A, prm);
  }

  /*!
   * \brief Solves the problem for the given right-hand side.
   *
   * \param rhs The right-hand side.
   * \param x   The solution.
   */
  template <class Vec1, class Vec2>
  void solve(const Vec1& rhs, Vec2& x) const
  {
    (*S)(rhs, x);
  }

 private:

  typedef DistributedDirectSolverBase<value_type, DistributedEigenSparseLUDirectSolver<value_type>> Base;
  params prm;
  std::shared_ptr<Solver> S;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::mpi::direct

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
