// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* BuiltinBackend.h                                            (C) 2026-2026 */
/*                                                                           */
/* Builtin backend using CSR matrix.                                         */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_BUILTINBACKEND_H
#define ARCANE_ALINA_BUILTINBACKEND_H
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

#pragma GCC diagnostic ignored "-Wconversion"

#include <vector>
#include <numeric>
#include <memory>
#include <random>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/alina/util.h>
#include <arcane/alina/ValueTypeInterface.h>
#include <arcane/alina/SkylineLUSolver.h>
#include <arcane/alina/DenseMatrixInverseImpl.h>
#include <arcane/alina/detail_sort_row.h>
#include <arcane/alina/SparseMatrixMatrixProduct.h>
#include <arcane/alina/backend_detail_matrix_ops.h>
#include <arcane/alina/CSRMatrix.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::backend
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Sort rows of the matrix column-wise.
template <typename V, typename C, typename P>
void sort_rows(CSRMatrix<V, C, P>& A)
{
  const size_t n = rows(A);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
    P beg = A.ptr[i];
    P end = A.ptr[i + 1];
    Alina::detail::sort_row(A.col + beg, A.val + beg, end - beg);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Transpose of a sparse matrix.
template <typename V, typename C, typename P>
std::shared_ptr<CSRMatrix<V, C, P>> transpose(const CSRMatrix<V, C, P>& A)
{
  const size_t n = rows(A);
  const size_t m = cols(A);
  const size_t nnz = nonzeros(A);

  auto T = std::make_shared<CSRMatrix<V, C, P>>();
  T->set_size(m, n, true);

  for (size_t j = 0; j < nnz; ++j)
    ++(T->ptr[A.col[j] + 1]);

  T->scan_row_sizes();
  T->set_nonzeros();

  for (size_t i = 0; i < n; i++) {
    for (P j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
      P head = T->ptr[A.col[j]]++;

      T->col[head] = static_cast<C>(i);
      T->val[head] = math::adjoint(A.val[j]);
    }
  }

  std::rotate(T->ptr, T->ptr + m, T->ptr + m + 1);
  T->ptr[0] = 0;

  return T;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Matrix-matrix product.
template <class Val, class Col, class Ptr>
std::shared_ptr<CSRMatrix<Val, Col, Ptr>>
product(const CSRMatrix<Val, Col, Ptr>& A, const CSRMatrix<Val, Col, Ptr>& B, bool sort = false)
{
  auto C = std::make_shared<CSRMatrix<Val, Col, Ptr>>();

#ifdef _OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

  if (nt > 16) {
    spgemm_rmerge(A, B, *C);
  }
  else {
    spgemm_saad(A, B, *C, sort);
  }

  return C;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Sum of two matrices
template <class Val, class Col, class Ptr>
std::shared_ptr<CSRMatrix<Val, Col, Ptr>>
sum(Val alpha, const CSRMatrix<Val, Col, Ptr>& A, Val beta,
    const CSRMatrix<Val, Col, Ptr>& B, bool sort = false)
{
  typedef ptrdiff_t Idx;

  auto C = std::make_shared<CSRMatrix<Val, Col, Ptr>>();
  precondition(A.nrows == B.nrows && A.ncols == B.ncols, "matrices should have same shape!");
  C->set_size(A.nrows, A.ncols);

  C->ptr[0] = 0;

#pragma omp parallel
  {
    std::vector<ptrdiff_t> marker(C->ncols, -1);

#pragma omp for
    for (Idx i = 0; i < static_cast<Idx>(C->nrows); ++i) {
      Idx C_cols = 0;

      for (Idx j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
        Idx c = A.col[j];

        if (marker[c] != i) {
          marker[c] = i;
          ++C_cols;
        }
      }

      for (Idx j = B.ptr[i], e = B.ptr[i + 1]; j < e; ++j) {
        Idx c = B.col[j];

        if (marker[c] != i) {
          marker[c] = i;
          ++C_cols;
        }
      }

      C->ptr[i + 1] = C_cols;
    }
  }

  C->set_nonzeros(C->scan_row_sizes());

#pragma omp parallel
  {
    std::vector<ptrdiff_t> marker(C->ncols, -1);

#pragma omp for
    for (Idx i = 0; i < static_cast<Idx>(C->nrows); ++i) {
      Idx row_beg = C->ptr[i];
      Idx row_end = row_beg;

      for (Idx j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
        Idx c = A.col[j];
        Val v = alpha * A.val[j];

        if (marker[c] < row_beg) {
          marker[c] = row_end;
          C->col[row_end] = c;
          C->val[row_end] = v;
          ++row_end;
        }
        else {
          C->val[marker[c]] += v;
        }
      }

      for (Idx j = B.ptr[i], e = B.ptr[i + 1]; j < e; ++j) {
        Idx c = B.col[j];
        Val v = beta * B.val[j];

        if (marker[c] < row_beg) {
          marker[c] = row_end;
          C->col[row_end] = c;
          C->val[row_end] = v;
          ++row_end;
        }
        else {
          C->val[marker[c]] += v;
        }
      }

      if (sort)
        Alina::detail::sort_row(C->col + row_beg, C->val + row_beg, row_end - row_beg);
    }
  }

  return C;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Scale matrix values.
template <class Val, class Col, class Ptr, class T>
void scale(CSRMatrix<Val, Col, Ptr>& A, T s)
{
  ptrdiff_t n = backend::rows(A);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n; ++i) {
    for (ptrdiff_t j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j)
      A.val[j] *= s;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Reduce matrix to a pointwise one
template <class value_type, class col_type, class ptr_type>
std::shared_ptr<CSRMatrix<typename math::scalar_of<value_type>::type, col_type, ptr_type>>
pointwise_matrix(const CSRMatrix<value_type, col_type, ptr_type>& A, unsigned block_size)
{
  typedef value_type V;
  typedef typename math::scalar_of<V>::type S;

  ARCANE_ALINA_TIC("pointwise_matrix");
  const ptrdiff_t n = A.nrows;
  const ptrdiff_t m = A.ncols;
  const ptrdiff_t np = n / block_size;
  const ptrdiff_t mp = m / block_size;

  precondition(np * block_size == n,
               "Matrix size should be divisible by block_size");

  auto ap = std::make_shared<CSRMatrix<S, col_type, ptr_type>>();
  auto& Ap = *ap;

  Ap.set_size(np, mp, true);

#pragma omp parallel
  {
    std::vector<ptr_type> j(block_size);
    std::vector<ptr_type> e(block_size);

    // Count number of nonzeros in block matrix.
#pragma omp for
    for (ptrdiff_t ip = 0; ip < np; ++ip) {
      ptrdiff_t ia = ip * block_size;
      col_type cur_col = 0;
      bool done = true;

      for (unsigned k = 0; k < block_size; ++k) {
        ptr_type beg = j[k] = A.ptr[ia + k];
        ptr_type end = e[k] = A.ptr[ia + k + 1];

        if (beg == end)
          continue;

        col_type c = A.col[beg];

        if (done) {
          done = false;
          cur_col = c;
        }
        else {
          cur_col = std::min(cur_col, c);
        }
      }

      while (!done) {
        cur_col /= block_size;
        ++Ap.ptr[ip + 1];

        done = true;
        col_type col_end = (cur_col + 1) * block_size;
        for (unsigned k = 0; k < block_size; ++k) {
          ptr_type beg = j[k];
          ptr_type end = e[k];

          while (beg < end) {
            col_type c = A.col[beg++];

            if (c >= col_end) {
              if (done) {
                done = false;
                cur_col = c;
              }
              else {
                cur_col = std::min(cur_col, c);
              }

              break;
            }
          }

          j[k] = beg;
        }
      }
    }
  }

  Ap.set_nonzeros(Ap.scan_row_sizes());

#pragma omp parallel
  {
    std::vector<ptr_type> j(block_size);
    std::vector<ptr_type> e(block_size);

#pragma omp for
    for (ptrdiff_t ip = 0; ip < np; ++ip) {
      ptrdiff_t ia = ip * block_size;
      col_type cur_col = 0;
      ptr_type head = Ap.ptr[ip];
      bool done = true;

      for (unsigned k = 0; k < block_size; ++k) {
        ptr_type beg = j[k] = A.ptr[ia + k];
        ptr_type end = e[k] = A.ptr[ia + k + 1];

        if (beg == end)
          continue;

        col_type c = A.col[beg];

        if (done) {
          done = false;
          cur_col = c;
        }
        else {
          cur_col = std::min(cur_col, c);
        }
      }

      while (!done) {
        cur_col /= block_size;

        Ap.col[head] = cur_col;

        done = true;
        bool first = true;
        S cur_val = math::zero<S>();

        col_type col_end = (cur_col + 1) * block_size;
        for (unsigned k = 0; k < block_size; ++k) {
          ptr_type beg = j[k];
          ptr_type end = e[k];

          while (beg < end) {
            col_type c = A.col[beg];
            S v = math::norm(A.val[beg]);
            ++beg;

            if (c >= col_end) {
              if (done) {
                done = false;
                cur_col = c;
              }
              else {
                cur_col = std::min(cur_col, c);
              }

              break;
            }

            if (first) {
              first = false;
              cur_val = v;
            }
            else {
              cur_val = std::max(cur_val, v);
            }
          }

          j[k] = beg;
        }

        Ap.val[head++] = cur_val;
      }
    }
  }

  ARCANE_ALINA_TOC("pointwise_matrix");
  return ap;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/** NUMA-aware vector container. */
template <class T>
class numa_vector
{
 public:

  typedef T value_type;

  numa_vector()
  : n(0)
  , p(0)
  {}

  numa_vector(size_t n, bool init = true)
  : n(n)
  , p(new T[n])
  {
    if (init) {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
        p[i] = math::zero<T>();
    }
  }

  void resize(size_t size, bool init = true)
  {
    delete[] p;
    p = 0;

    n = size;
    p = new T[n];

    if (init) {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
        p[i] = math::zero<T>();
    }
  }

  template <class Vector>
  numa_vector(const Vector& other,
              typename std::enable_if<!std::is_integral<Vector>::value, int>::type = 0)
  : n(other.size())
  , p(new T[n])
  {
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
      p[i] = other[i];
  }

  template <class Iterator>
  numa_vector(Iterator beg, Iterator end)
  : n(std::distance(beg, end))
  , p(new T[n])
  {
    static_assert(std::is_same<
                  std::random_access_iterator_tag,
                  typename std::iterator_traits<Iterator>::iterator_category>::value,
                  "Iterator has to be random access");

#pragma omp parallel for
    for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
      p[i] = beg[i];
  }

  ~numa_vector()
  {
    delete[] p;
    p = 0;
  }

  inline size_t size() const
  {
    return n;
  }

  inline const T& operator[](size_t i) const
  {
    return p[i];
  }

  inline T& operator[](size_t i)
  {
    return p[i];
  }

  inline const T* data() const
  {
    return p;
  }

  inline T* data()
  {
    return p;
  }

  void swap(numa_vector& other)
  {
    std::swap(n, other.n);
    std::swap(p, other.p);
  }

 private:

  size_t n;
  T* p;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Diagonal of a matrix
template <typename V, typename C, typename P>
std::shared_ptr<numa_vector<V>> diagonal(const CSRMatrix<V, C, P>& A, bool invert = false)
{
  const size_t n = rows(A);
  auto dia = std::make_shared<numa_vector<V>>(n, false);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
    for (auto a = A.row_begin(i); a; ++a) {
      if (a.col() == i) {
        V d = a.value();
        if (invert) {
          d = math::is_zero(d) ? math::identity<V>() : math::inverse(d);
        }
        (*dia)[i] = d;
        break;
      }
    }
  }

  return dia;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.
template <bool scale, class Matrix>
static typename math::scalar_of<typename backend::value_type<Matrix>::type>::type
spectral_radius(const Matrix& A, int power_iters = 0)
{
  ARCANE_ALINA_TIC("spectral radius");
  typedef typename backend::value_type<Matrix>::type value_type;
  typedef typename math::rhs_of<value_type>::type rhs_type;
  typedef typename math::scalar_of<value_type>::type scalar_type;

  const ptrdiff_t n = backend::rows(A);
  scalar_type radius;

  if (power_iters <= 0) {
    // Use Gershgorin disk theorem.
    radius = 0;

#pragma omp parallel
    {
      scalar_type emax = 0;
      value_type dia = math::identity<value_type>();

#pragma omp for nowait
      for (ptrdiff_t i = 0; i < n; ++i) {
        scalar_type s = 0;

        for (ptrdiff_t j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
          ptrdiff_t c = A.col[j];
          value_type v = A.val[j];

          s += math::norm(v);

          if (scale && c == i)
            dia = v;
        }

        if (scale)
          s *= math::norm(math::inverse(dia));

        emax = std::max(emax, s);
      }

#pragma omp critical
      radius = std::max(radius, emax);
    }
  }
  else {
    // Power method.
    backend::numa_vector<rhs_type> b0(n, false), b1(n, false);

    // Fill the initial vector with random values.
    // Also extract the inverted matrix diagonal values.
    scalar_type b0_norm = 0;
#pragma omp parallel
    {
#ifdef _OPENMP
      int tid = omp_get_thread_num();
#else
      int tid = 0;
#endif
      std::mt19937 rng(tid);
      std::uniform_real_distribution<scalar_type> rnd(-1, 1);

      scalar_type loc_norm = 0;

#pragma omp for nowait
      for (ptrdiff_t i = 0; i < n; ++i) {
        rhs_type v = math::constant<rhs_type>(rnd(rng));

        b0[i] = v;
        loc_norm += math::norm(math::inner_product(v, v));
      }

#pragma omp critical
      b0_norm += loc_norm;
    }

    // Normalize b0
    b0_norm = 1 / sqrt(b0_norm);
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n; ++i) {
      b0[i] = b0_norm * b0[i];
    }

    for (int iter = 0; iter < power_iters;) {
      // b1 = scale ? (D^1 * A) * b0 : A * b0
      // b1_norm = ||b1||
      // radius = <b1,b0>
      scalar_type b1_norm = 0;
      radius = 0;
#pragma omp parallel
      {
        scalar_type loc_norm = 0;
        scalar_type loc_radi = 0;
        value_type dia = math::identity<value_type>();

#pragma omp for nowait
        for (ptrdiff_t i = 0; i < n; ++i) {
          rhs_type s = math::zero<rhs_type>();

          for (ptrdiff_t j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
            ptrdiff_t c = A.col[j];
            value_type v = A.val[j];
            if (scale && c == i)
              dia = v;
            s += v * b0[c];
          }

          if (scale)
            s = math::inverse(dia) * s;

          loc_norm += math::norm(math::inner_product(s, s));
          loc_radi += math::norm(math::inner_product(s, b0[i]));

          b1[i] = s;
        }

#pragma omp critical
        {
          b1_norm += loc_norm;
          radius += loc_radi;
        }
      }

      if (++iter < power_iters) {
        // b0 = b1 / b1_norm
        b1_norm = 1 / sqrt(b1_norm);
#pragma omp parallel for
        for (ptrdiff_t i = 0; i < n; ++i) {
          b0[i] = b1_norm * b1[i];
        }
      }
    }
  }
  ARCANE_ALINA_TOC("spectral radius");

  return radius < 0 ? static_cast<scalar_type>(2) : radius;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/**
 * The builtin backend does not have any dependencies, and uses OpenMP for
 * parallelization. Matrices are stored in the CRS format, and vectors are
 * instances of ``std::vector<value_type>``. There is no usual overhead of
 * moving the constructed hierarchy to the builtin backend, since the backend
 * is used internally during setup.
 */
template <typename ValueType, typename ColumnType = ptrdiff_t, typename PointerType = ColumnType>
struct BuiltinBackend
{
  typedef ValueType value_type;
  typedef ColumnType index_type;
  typedef ColumnType col_type;
  typedef PointerType ptr_type;

  typedef typename math::rhs_of<value_type>::type rhs_type;

  struct provides_row_iterator : std::true_type
  {};

  typedef CSRMatrix<value_type, col_type, ptr_type> matrix;
  typedef numa_vector<rhs_type> vector;
  typedef numa_vector<value_type> matrix_diagonal;
  typedef solver::SkylineLUSolver<value_type> direct_solver;

  /// The backend has no parameters.
  typedef Alina::detail::empty_params params;

  static std::string name() { return "builtin"; }

  // Copy matrix. This is a noop for builtin backend.
  static std::shared_ptr<matrix>
  copy_matrix(std::shared_ptr<matrix> A, const params&)
  {
    return A;
  }

  // Copy vector to builtin backend.
  template <class T>
  static std::shared_ptr<numa_vector<T>>
  copy_vector(const std::vector<T>& x, const params&)
  {
    return std::make_shared<numa_vector<T>>(x);
  }

  // Copy vector to builtin backend. This is a noop for builtin backend.
  template <class T>
  static std::shared_ptr<numa_vector<T>>
  copy_vector(std::shared_ptr<numa_vector<T>> x, const params&)
  {
    return x;
  }

  // Create vector of the specified size.
  static std::shared_ptr<vector>
  create_vector(size_t size, const params&)
  {
    return std::make_shared<vector>(size);
  }

  struct gather
  {
    std::vector<ptrdiff_t> I;

    gather(size_t /*size*/, const std::vector<ptrdiff_t>& I, const params&)
    : I(I)
    {}

    template <class InVec, class OutVec>
    void operator()(const InVec& vec, OutVec& vals) const
    {
      for (size_t i = 0; i < I.size(); ++i)
        vals[i] = vec[I[i]];
    }
  };

  struct scatter
  {
    std::vector<ptrdiff_t> I;

    scatter(size_t /*size*/, const std::vector<ptrdiff_t>& I, const params&)
    : I(I)
    {}

    template <class InVec, class OutVec>
    void operator()(const InVec& vals, OutVec& vec) const
    {
      for (size_t i = 0; i < I.size(); ++i)
        vec[I[i]] = vals[i];
    }
  };

  // Create direct solver for coarse level
  static std::shared_ptr<direct_solver>
  create_solver(std::shared_ptr<matrix> A, const params&)
  {
    return std::make_shared<direct_solver>(*A);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class T>
struct is_builtin_vector : std::false_type
{};

template <class V>
struct is_builtin_vector<std::vector<V>> : std::is_arithmetic<V>
{};

template <class V>
struct is_builtin_vector<numa_vector<V>> : std::true_type
{};

//---------------------------------------------------------------------------
// Specialization of backend interface
//---------------------------------------------------------------------------
template <typename T1, typename T2>
struct backends_compatible<BuiltinBackend<T1>, BuiltinBackend<T2>> : std::true_type
{};

template <typename V, typename C, typename P>
struct rows_impl<CSRMatrix<V, C, P>>
{
  static size_t get(const CSRMatrix<V, C, P>& A)
  {
    return A.nrows;
  }
};

template <typename V, typename C, typename P>
struct cols_impl<CSRMatrix<V, C, P>>
{
  static size_t get(const CSRMatrix<V, C, P>& A)
  {
    return A.ncols;
  }
};

template <class Vec>
struct bytes_impl<Vec, typename std::enable_if<is_builtin_vector<Vec>::value>::type>
{
  static size_t get(const Vec& x)
  {
    typedef typename backend::value_type<Vec>::type V;
    return x.size() * sizeof(V);
  }
};

template <typename V, typename C, typename P>
struct ptr_data_impl<CSRMatrix<V, C, P>>
{
  typedef const P* type;
  static type get(const CSRMatrix<V, C, P>& A)
  {
    return &A.ptr[0];
  }
};

template <typename V, typename C, typename P>
struct col_data_impl<CSRMatrix<V, C, P>>
{
  typedef const C* type;
  static type get(const CSRMatrix<V, C, P>& A)
  {
    return &A.col[0];
  }
};

template <typename V, typename C, typename P>
struct val_data_impl<CSRMatrix<V, C, P>>
{
  typedef const V* type;
  static type get(const CSRMatrix<V, C, P>& A)
  {
    return &A.val[0];
  }
};

template <typename V, typename C, typename P>
struct nonzeros_impl<CSRMatrix<V, C, P>>
{
  static size_t get(const CSRMatrix<V, C, P>& A)
  {
    return A.nrows == 0 ? 0 : A.ptr[A.nrows];
  }
};

template <typename V, typename C, typename P>
struct row_nonzeros_impl<CSRMatrix<V, C, P>>
{
  static size_t get(const CSRMatrix<V, C, P>& A, size_t row)
  {
    return A.ptr[row + 1] - A.ptr[row];
  }
};

template <class Vec>
struct clear_impl<Vec, typename std::enable_if<is_builtin_vector<Vec>::value>::type>
{
  static void apply(Vec& x)
  {
    typedef typename backend::value_type<Vec>::type V;

    const size_t n = x.size();
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
      x[i] = math::zero<V>();
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Vec1, class Vec2>
struct inner_product_impl<Vec1, Vec2,
                          typename std::enable_if<
                          is_builtin_vector<Vec1>::value &&
                          is_builtin_vector<Vec2>::value>::type>
{
  typedef typename value_type<Vec1>::type V;

  typedef typename math::inner_product_impl<V>::return_type return_type;

  static return_type get(const Vec1& x, const Vec2& y)
  {
#ifdef _OPENMP
    if (omp_get_max_threads() > 1) {
      return parallel(x, y);
    }
    else
#endif
    {
      return serial(x, y);
    }
  }

  static return_type serial(const Vec1& x, const Vec2& y)
  {
    const size_t n = x.size();

    return_type s = math::zero<return_type>();
    return_type c = math::zero<return_type>();

    for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
      return_type d = math::inner_product(x[i], y[i]) - c;
      return_type t = s + d;
      c = (t - s) - d;
      s = t;
    }

    return s;
  }

#ifdef _OPENMP
#ifndef ARCANE_ALINA_MAX_OPENMP_THREADS
#define ARCANE_ALINA_MAX_OPENMP_THREADS 64
#endif
  static return_type parallel(const Vec1& x, const Vec2& y)
  {
    const size_t n = x.size();
    return_type _sum_stat[ARCANE_ALINA_MAX_OPENMP_THREADS];
    std::vector<return_type> _sum_dyna;
    return_type* sum;

    const int nt = omp_get_max_threads();

    if (nt < 64) {
      sum = _sum_stat;
      for (int i = 0; i < nt; ++i) {
        sum[i] = math::zero<return_type>();
      }
    }
    else {
      _sum_dyna.resize(nt, math::zero<return_type>());
      sum = _sum_dyna.data();
    }

#pragma omp parallel
    {
      const int tid = omp_get_thread_num();

      return_type s = math::zero<return_type>();
      return_type c = math::zero<return_type>();

#pragma omp for nowait
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        return_type d = math::inner_product(x[i], y[i]) - c;
        return_type t = s + d;
        c = (t - s) - d;
        s = t;
      }

      sum[tid] = s;
    }

    return std::accumulate(sum, sum + nt, math::zero<return_type>());
  }
#endif
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class A, class Vec1, class B, class Vec2>
struct axpby_impl<A, Vec1, B, Vec2, typename std::enable_if<is_builtin_vector<Vec1>::value && is_builtin_vector<Vec2>::value>::type>
{
  static void apply(A a, const Vec1& x, B b, Vec2& y)
  {
    const size_t n = x.size();
    if (!math::is_zero(b)) {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        y[i] = a * x[i] + b * y[i];
      }
    }
    else {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        y[i] = a * x[i];
      }
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class A, class Vec1, class B, class Vec2, class C, class Vec3>
struct axpbypcz_impl<A, Vec1, B, Vec2, C, Vec3,
                     typename std::enable_if<
                     is_builtin_vector<Vec1>::value &&
                     is_builtin_vector<Vec2>::value &&
                     is_builtin_vector<Vec3>::value>::type>
{
  static void apply(A a, const Vec1& x, B b, const Vec2& y, C c, Vec3& z)
  {
    const size_t n = x.size();
    if (!math::is_zero(c)) {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        z[i] = a * x[i] + b * y[i] + c * z[i];
      }
    }
    else {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        z[i] = a * x[i] + b * y[i];
      }
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Alpha, class Vec1, class Vec2, class Beta, class Vec3>
struct vmul_impl<Alpha, Vec1, Vec2, Beta, Vec3,
                 typename std::enable_if<
                 is_builtin_vector<Vec1>::value &&
                 is_builtin_vector<Vec2>::value &&
                 is_builtin_vector<Vec3>::value &&
                 math::static_rows<typename value_type<Vec1>::type>::value == math::static_rows<typename value_type<Vec2>::type>::value &&
                 math::static_rows<typename value_type<Vec1>::type>::value == math::static_rows<typename value_type<Vec3>::type>::value>::type>
{
  static void apply(Alpha a, const Vec1& x, const Vec2& y, Beta b, Vec3& z)
  {
    const size_t n = x.size();
    if (!math::is_zero(b)) {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        z[i] = a * x[i] * y[i] + b * z[i];
      }
    }
    else {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        z[i] = a * x[i] * y[i];
      }
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Support for mixed scalar/nonscalar types
template <class Alpha, class Vec1, class Vec2, class Beta, class Vec3>
struct vmul_impl<Alpha, Vec1, Vec2, Beta, Vec3,
                 typename std::enable_if<is_builtin_vector<Vec1>::value &&
                                         is_builtin_vector<Vec2>::value &&
                                         is_builtin_vector<Vec3>::value &&
                                         (math::static_rows<typename value_type<Vec1>::type>::value != math::static_rows<typename value_type<Vec2>::type>::value ||
                                          math::static_rows<typename value_type<Vec1>::type>::value != math::static_rows<typename value_type<Vec3>::type>::value)>::type>
{
  static void apply(Alpha a, const Vec1& x, const Vec2& y, Beta b, Vec3& z)
  {
    typedef typename value_type<Vec1>::type M_type;
    auto Y = backend::reinterpret_as_rhs<M_type>(y);
    auto Z = backend::reinterpret_as_rhs<M_type>(z);

    const size_t n = x.size();

    if (!math::is_zero(b)) {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        Z[i] = a * x[i] * Y[i] + b * Z[i];
      }
    }
    else {
#pragma omp parallel for
      for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        Z[i] = a * x[i] * Y[i];
      }
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Vec1, class Vec2>
struct copy_impl<Vec1, Vec2,
                 typename std::enable_if<
                 is_builtin_vector<Vec1>::value &&
                 is_builtin_vector<Vec2>::value>::type>
{
  static void apply(const Vec1& x, Vec2& y)
  {
    const size_t n = x.size();
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
      y[i] = x[i];
    }
  }
};

template <class MatrixValue, class Vector, bool IsConst>
struct reinterpret_as_rhs_impl<MatrixValue, Vector, IsConst,
typename std::enable_if<is_builtin_vector<Vector>::value>::type>
{
  typedef typename backend::value_type<Vector>::type src_type;
  typedef typename math::scalar_of<src_type>::type scalar_type;
  typedef typename math::rhs_of<MatrixValue>::type rhs_type;
  typedef typename math::replace_scalar<rhs_type, scalar_type>::type dst_type;
  typedef typename std::conditional<IsConst, const dst_type*, dst_type*>::type ptr_type;
  typedef iterator_range<ptr_type> return_type;

  template <class V>
  static return_type get(V&& x)
  {
    auto ptr = reinterpret_cast<ptr_type>(&x[0]);
    const size_t n = x.size() * sizeof(src_type) / sizeof(dst_type);
    return make_iterator_range(ptr, ptr + n);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace detail
{

template <typename V, typename C, typename P>
struct use_builtin_matrix_ops<CSRMatrix<V, C, P>>
: std::true_type
{};

} // namespace detail

} // namespace Arcane::Alina::backend

// Allow to use boost::iterator_range as vector in builtin backend:
namespace boost
{
  template <class Iterator> class iterator_range;
}

namespace Arcane::Alina::backend
{
  template <class Iterator>
  struct is_builtin_vector< Alina::iterator_range<Iterator> > : std::true_type {};

  template <class Iterator>
  struct is_builtin_vector< boost::iterator_range<Iterator> > : std::true_type {};
}

namespace Arcane::Alina::detail
{

// Backend with scalar value_type of highest precision.
template <class V1, class V2>
struct common_scalar_backend<backend::BuiltinBackend<V1>, backend::BuiltinBackend<V2>,
                             typename std::enable_if<math::static_rows<V1>::value != 1 || math::static_rows<V2>::value != 1>::type>
{
  typedef typename math::scalar_of<V1>::type S1;
  typedef typename math::scalar_of<V2>::type S2;

  typedef typename std::conditional<(sizeof(S1) > sizeof(S2)), backend::BuiltinBackend<S1>, backend::BuiltinBackend<S2>>::type type;
};

} // namespace detail

#endif
