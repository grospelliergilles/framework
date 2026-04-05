// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MPRelaxation.h                                              (C) 2000-2026 */
/*                                                                           */
/* Relaxtion with message passing support.                                   */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MP_RELAXATION_H
#define ARCANE_ALINA_MP_RELAXATION_H
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

#include <arcane/alina/relaxation.h>

#include <arcane/alina/backend_interface.h>
#include <arcane/alina/backend_builtin.h>

#include <arcane/alina/mpi/mp_distributed_matrix.h>
#include <arcane/alina/mpi/mp_util.h>

#include <memory>
#include <vector>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::mpi::relaxation
{

template <class Backend>
struct chebyshev
: public Alina::relaxation::chebyshev<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::chebyshev<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  chebyshev(const distributed_matrix<Backend>& A,
            const params& prm = params(),
            const backend_params& bprm = backend_params())
  : Base(A, prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct damped_jacobi
: public Alina::relaxation::damped_jacobi<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::damped_jacobi<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  damped_jacobi(const distributed_matrix<Backend>& A,
                const params& prm = params(),
                const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct gauss_seidel
: public Alina::relaxation::gauss_seidel<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::gauss_seidel<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  gauss_seidel(const distributed_matrix<Backend>& A,
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
struct ilu0
: public Alina::relaxation::ilu0<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ilu0<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  ilu0(const distributed_matrix<Backend>& A,
       const params& prm = params(),
       const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct iluk : public Alina::relaxation::iluk<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::iluk<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  iluk(const distributed_matrix<Backend>& A,
       const params& prm = params(),
       const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct ilup : public Alina::relaxation::ilup<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ilup<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  ilup(const distributed_matrix<Backend>& A,
       const params& prm = params(),
       const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct ilut : public Alina::relaxation::ilut<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::ilut<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  ilut(const distributed_matrix<Backend>& A,
       const params& prm = params(),
       const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct spai0
{
  typedef Backend backend_type;
  typedef typename Backend::value_type value_type;
  typedef typename Backend::matrix_diagonal matrix_diagonal;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  typedef Alina::detail::empty_params params;
  typedef typename Backend::params backend_params;

  spai0(const distributed_matrix<Backend>& A,
        const params&, const backend_params& bprm = backend_params())
  {
    typedef backend::crs<value_type> build_matrix;

    const ptrdiff_t n = A.loc_rows();
    const build_matrix& A_loc = *A.local();
    const build_matrix& A_rem = *A.remote();

    auto m = std::make_shared<backend::numa_vector<value_type>>(n, false);
    typedef backend::crs<value_type> build_matrix;

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
  void apply_pre(
  const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
  {
    static const scalar_type one = math::identity<scalar_type>();
    backend::residual(rhs, A, x, tmp);
    backend::vmul(one, *M, tmp, one, x);
  }

  /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
  template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
  void apply_post(
  const Matrix& A, const VectorRHS& rhs, VectorX& x, VectorTMP& tmp) const
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
struct spai1 : public Alina::relaxation::spai1<Backend>
{
  typedef Backend backend_type;
  typedef Alina::relaxation::spai1<Backend> Base;
  typedef typename Backend::params backend_params;
  typedef typename Base::params params;

  spai1(
  const distributed_matrix<Backend>& A,
  const params& prm = params(),
  const backend_params& bprm = backend_params())
  : Base(*A.local(), prm, bprm)
  {}
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Relaxation>
struct as_preconditioner
{
  typedef typename Relaxation::params params;
  typedef typename Relaxation::backend_type backend_type;
  typedef typename backend_type::params backend_params;
  typedef typename backend_type::value_type value_type;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  typedef distributed_matrix<backend_type> matrix;
  typedef typename backend_type::vector vector;

  template <class Matrix>
  as_preconditioner(communicator comm,
                    const Matrix& A,
                    const params& prm = params(),
                    const backend_params& bprm = backend_params())
  : A(std::make_shared<matrix>(comm, A, backend::rows(A)))
  , S(A, prm, bprm)
  {
    this->A->move_to_backend(bprm);
  }

  as_preconditioner(communicator,
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

  friend std::ostream& operator<<(std::ostream& os, const as_preconditioner& p)
  {
    os << "Relaxation as preconditioner" << std::endl;
    os << "  unknowns: " << p.system_matrix().glob_rows() << std::endl;
    os << "  nonzeros: " << p.system_matrix().glob_nonzeros() << std::endl;

    return os;
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace amgcl::mpi::relaxation

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
