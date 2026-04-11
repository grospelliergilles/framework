// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedRelaxation.h                                     (C) 2000-2026 */
/*                                                                           */
/* Relaxation with distribution support.                                     */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_DISTRIBUTEDRELAXATION_H
#define ARCANE_ALINA_DISTRIBUTEDRELAXATION_H
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
#include <arcane/alina/BuiltinBackend.h>
#include <arcane/alina/Relaxation.h>
#include <arcane/alina/DistributedMatrix.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedChebyshevRelaxation
: public Alina::relaxation::ChebyshevRelaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ChebyshevRelaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedChebyshevRelaxation(const DistributedMatrix<Backend>& A,
                                 const params& prm = params(),
                                 const backend_params& bprm = backend_params())
  : Base(A, prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedDampedJacobiRelaxation
: public Alina::relaxation::DampedJacobiRelaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::DampedJacobiRelaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedDampedJacobiRelaxation(const DistributedMatrix<Backend>& A,
                                    const params& prm = params(),
                                    const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedGaussSeidelRelaxation
: public Alina::relaxation::GaussSeidelRelaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::GaussSeidelRelaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedGaussSeidelRelaxation(const DistributedMatrix<Backend>& A,
                                   const params& prm = params(),
                                   const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}

  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_pre(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& t) const
  {
    Base::apply_pre(*A.local_backend(), rhs, x, t);
  }

  /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_post(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& t) const
  {
    Base::apply_post(*A.local_backend(), rhs, x, t);
  }

  template <class Matrix, class VectorRHS, class VectorX>
  void apply(const Matrix& A, const VectorRHS& rhs, VectorX& x) const
  {
    Base::apply(*A.local_backend(), rhs, x);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedILU0Relaxation
: public Alina::relaxation::ILU0Relaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ILU0Relaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedILU0Relaxation(const DistributedMatrix<Backend>& A,
                            const params& prm = params(),
                            const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedILUKRelaxation
: public Alina::relaxation::ILUKRelaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ILUKRelaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedILUKRelaxation(const DistributedMatrix<Backend>& A,
                            const params& prm = params(),
                            const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedILUPRelaxation
: public Alina::relaxation::ILUPRelaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ILUPRelaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedILUPRelaxation(const DistributedMatrix<Backend>& A,
                            const params& prm = params(),
                            const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedILUTRelaxation : public Alina::relaxation::ILUTRelaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ILUTRelaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedILUTRelaxation(const DistributedMatrix<Backend>& A,
                            const params& prm = params(),
                            const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedSPAI0Relaxation
{
  typedef Backend backend_type;
  typedef typename Backend::value_type value_type;
  typedef typename Backend::matrix_diagonal matrix_diagonal;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  typedef Alina::detail::empty_params params;
  typedef typename Backend::params backend_params;

  DistributedSPAI0Relaxation(const DistributedMatrix<Backend>& A,
                             const params&, const backend_params& bprm = backend_params())
  {
    typedef backend::CSRMatrix<value_type> build_matrix;

    const ptrdiff_t n = A.loc_rows();
    const build_matrix& A_loc = *A.local();
    const build_matrix& A_rem = *A.remote();

    auto m = std::make_shared<backend::numa_vector<value_type>>(n, false);
    typedef backend::CSRMatrix<value_type> build_matrix;

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n; ++i) {
      value_type num = math::zero<value_type>();
      scalar_type den = math::zero<scalar_type>();

      for (ptrdiff_t j = A_loc.ptr[i], e = A_loc.ptr[i + 1]; j < e; ++j) {
        value_type v = A_loc.val[j];
        scalar_type norm_v = math::norm(v);
        den += norm_v * norm_v;
        if (A_loc.col[j] == i)
          num += v;
      }

      for (ptrdiff_t j = A_rem.ptr[i], e = A_rem.ptr[i + 1]; j < e; ++j) {
        value_type v = A_rem.val[j];
        scalar_type norm_v = math::norm(v);
        den += norm_v * norm_v;
      }

      (*m)[i] = math::inverse(den) * num;
    }

    M = Backend::copy_vector(m, bprm);
  }

  /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_pre(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    static const scalar_type one = math::identity<scalar_type>();
    backend::residual(rhs, A, x, tmp);
    backend::vmul(one, *M, tmp, one, x);
  }

  /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_post(const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    static const scalar_type one = math::identity<scalar_type>();
    backend::residual(rhs, A, x, tmp);
    backend::vmul(one, *M, tmp, one, x);
  }

  template <class Matrix, class VectorRHS, class VectorX>
  void apply(const Matrix&, const VectorRHS& rhs, VectorX& x) const
  {
    backend::vmul(math::identity<scalar_type>(), *M, rhs, math::zero<scalar_type>(), x);
  }

 private:

  std::shared_ptr<matrix_diagonal> M;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedSPAI1Relaxation
: public Alina::relaxation::SPAI1Relaxation<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::SPAI1Relaxation<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  DistributedSPAI1Relaxation(const DistributedMatrix<Backend>& A,
                             const params& prm = params(),
                             const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Use a relaxation as a distributed preconditioner.
 */
template <class Relaxation>
struct AsDistributedPreconditioner
{
  typedef typename Relaxation::params params;
  typedef typename Relaxation::BackendType backend_type;
  using BackendType = backend_type;
  typedef typename backend_type::params backend_params;
  typedef typename backend_type::value_type value_type;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  typedef DistributedMatrix<backend_type> matrix;
  typedef typename backend_type::vector vector;

  template <class Matrix>
  AsDistributedPreconditioner(mpi_communicator comm,
                              const Matrix& A,
                              const params& prm = params(),
                              const backend_params& bprm = backend_params())
  : A(std::make_shared<matrix>(comm, A, backend::rows(A)))
  , S(A, prm, bprm)
  {
    this->A->move_to_backend(bprm);
  }

  AsDistributedPreconditioner(mpi_communicator,
                              std::shared_ptr<matrix> A,
                              const params& prm = params(),
                              const backend_params& bprm = backend_params())
  : A(A)
  , S(*A, prm, bprm)
  {
    this->A->move_to_backend(bprm);
  }

  template <class Vec1, class Vec2>
  void apply(const Vec1& rhs, Vec2&& x) const
  {
    S.apply(*A, rhs, x);
  }

  std::shared_ptr<matrix> system_matrix_ptr() const
  {
    return A;
  }

  const matrix& system_matrix() const
  {
    return *system_matrix_ptr();
  }

 private:

  std::shared_ptr<matrix> A;
  Relaxation S;

  friend std::ostream& operator<<(std::ostream& os, const AsDistributedPreconditioner& p)
  {
    os << "Relaxation as preconditioner" << std::endl;
    os << "  unknowns: " << p.system_matrix().glob_rows() << std::endl;
    os << "  nonzeros: " << p.system_matrix().glob_nonzeros() << std::endl;

    return os;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
