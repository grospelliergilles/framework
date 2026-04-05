// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Adapters.h                                                  (C) 2026-2026 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_ADAPTERS_H
#define ARCANE_ALINA_ADAPTERS_H
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

#include <type_traits>
#include <vector>
#include <numeric>
#include <tuple>

#include <boost/range/iterator_range.hpp>
#include <boost/range/size.hpp>
#include <boost/iterator/permutation_iterator.hpp>

#include <arcane/alina/util.h>
#include <arcane/alina/backend_builtin.h>
#include <arcane/alina/value_type_backend_interface.h>
#include <arcane/alina/backend_detail_matrix_ops.h>
#include <arcane/alina/reorder_cuthill_mckee.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::backend
{

//---------------------------------------------------------------------------
// Specialization of matrix interface
//---------------------------------------------------------------------------
template <typename N, typename PRng, typename CRng, typename VRng>
struct value_type<std::tuple<N, PRng, CRng, VRng>>
{
  typedef typename std::decay<decltype(std::declval<VRng>()[0])>::type type;
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct rows_impl<std::tuple<N, PRng, CRng, VRng>>
{
  static size_t get(const std::tuple<N, PRng, CRng, VRng>& A)
  {
    return std::get<0>(A);
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct cols_impl<std::tuple<N, PRng, CRng, VRng>>
{
  static size_t get(const std::tuple<N, PRng, CRng, VRng>& A)
  {
    return std::get<0>(A);
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct nonzeros_impl<std::tuple<N, PRng, CRng, VRng>>
{
  static size_t get(const std::tuple<N, PRng, CRng, VRng>& A)
  {
    return std::get<1>(A)[std::get<0>(A)];
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct row_iterator<std::tuple<N, PRng, CRng, VRng>>
{
  class type
  {
   public:

    typedef typename std::decay<decltype(std::declval<CRng>()[0])>::type col_type;
    typedef typename std::decay<decltype(std::declval<VRng>()[0])>::type val_type;

    type(const std::tuple<N, PRng, CRng, VRng>& A, size_t row)
    : m_col(std::begin(std::get<2>(A)))
    , m_end(std::begin(std::get<2>(A)))
    , m_val(std::begin(std::get<3>(A)))
    {
      typedef typename std::decay<decltype(std::declval<PRng>()[0])>::type ptr_type;

      ptr_type row_begin = std::get<1>(A)[row];
      ptr_type row_end = std::get<1>(A)[row + 1];

      m_col += row_begin;
      m_end += row_end;
      m_val += row_begin;
    }

    operator bool() const
    {
      return m_col != m_end;
    }

    type& operator++()
    {
      ++m_col;
      ++m_val;
      return *this;
    }

    col_type col() const
    {
      return *m_col;
    }

    val_type value() const
    {
      return *m_val;
    }

   private:

    typedef decltype(std::begin(std::declval<VRng>())) val_iterator;
    typedef decltype(std::begin(std::declval<CRng>())) col_iterator;

    col_iterator m_col;
    col_iterator m_end;
    val_iterator m_val;
  };
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct row_begin_impl<std::tuple<N, PRng, CRng, VRng>>
{
  typedef std::tuple<N, PRng, CRng, VRng> Matrix;
  static typename row_iterator<Matrix>::type
  get(const Matrix& matrix, size_t row)
  {
    return typename row_iterator<Matrix>::type(matrix, row);
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct row_nonzeros_impl<std::tuple<N, PRng, CRng, VRng>>
{
  typedef std::tuple<N, PRng, CRng, VRng> Matrix;

  static size_t get(const Matrix& A, size_t row)
  {
    return std::get<1>(A)[row + 1] - std::get<1>(A)[row];
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct ptr_data_impl<std::tuple<N, PRng, CRng, VRng>>
{
  typedef std::tuple<N, PRng, CRng, VRng> Matrix;
  typedef typename std::decay<decltype(std::declval<PRng>()[0])>::type ptr_type;
  typedef const ptr_type* type;
  static type get(const Matrix& A)
  {
    return &std::get<1>(A)[0];
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct col_data_impl<std::tuple<N, PRng, CRng, VRng>>
{
  typedef std::tuple<N, PRng, CRng, VRng> Matrix;
  typedef typename std::decay<decltype(std::declval<CRng>()[0])>::type col_type;
  typedef const col_type* type;
  static type get(const Matrix& A)
  {
    return &std::get<2>(A)[0];
  }
};

template <typename N, typename PRng, typename CRng, typename VRng>
struct val_data_impl<std::tuple<N, PRng, CRng, VRng>>
{
  typedef std::tuple<N, PRng, CRng, VRng> Matrix;
  typedef typename std::decay<decltype(std::declval<VRng>()[0])>::type val_type;
  typedef const val_type* type;
  static type get(const Matrix& A)
  {
    return &std::get<3>(A)[0];
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::backend

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::adapter
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Matrix, class BlockType>
struct block_matrix_adapter
{
  typedef BlockType value_type;
  static const int BlockSize = math::static_rows<BlockType>::value;

  const Matrix& A;

  block_matrix_adapter(const Matrix& A)
  : A(A)
  {
    precondition(
    backend::rows(A) % BlockSize == 0 &&
    backend::cols(A) % BlockSize == 0,
    "Matrix size is not divisible by block size!");
  }

  size_t rows() const
  {
    return backend::rows(A) / BlockSize;
  }

  size_t cols() const
  {
    return backend::cols(A) / BlockSize;
  }

  size_t nonzeros() const
  {
    // Just an estimate:
    return backend::nonzeros(A) / (BlockSize * BlockSize);
  }

  struct row_iterator
  {
    typedef typename backend::row_iterator<Matrix>::type Base;
    typedef ptrdiff_t col_type;
    typedef BlockType val_type;

    std::array<char, sizeof(Base) * BlockSize> buf;
    Base* base;

    bool done;
    col_type cur_col;
    val_type cur_val;

    row_iterator(const Matrix& A, col_type row)
    : done(true)
    {
      base = reinterpret_cast<Base*>(buf.data());
      for (int i = 0; i < BlockSize; ++i) {
        new (base + i) Base(backend::row_begin(A, row * BlockSize + i));

        if (base[i]) {
          col_type col = base[i].col() / BlockSize;
          if (done) {
            cur_col = col;
            done = false;
          }
          else {
            cur_col = std::min<col_type>(cur_col, col);
          }
        }
      }

      if (done)
        return;

      // While we are gathering the current value,
      // base iteratirs are advanced to the next block-column.
      cur_val = math::zero<val_type>();
      col_type end = (cur_col + 1) * BlockSize;
      for (int i = 0; i < BlockSize; ++i) {
        for (; base[i] && static_cast<ptrdiff_t>(base[i].col()) < end; ++base[i]) {
          cur_val(i, base[i].col() % BlockSize) = base[i].value();
        }
      }
    }

    ~row_iterator()
    {
      for (int i = 0; i < BlockSize; ++i)
        base[i].~Base();
    }

    operator bool() const
    {
      return !done;
    }

    row_iterator& operator++()
    {
      // Base iterators are already at the next block-column.
      // We just need to gather the current column and value.
      done = true;

      col_type end = (cur_col + 1) * BlockSize;
      for (int i = 0; i < BlockSize; ++i) {
        if (base[i]) {
          col_type col = base[i].col() / BlockSize;
          if (done) {
            cur_col = col;
            done = false;
          }
          else {
            cur_col = std::min<col_type>(cur_col, col);
          }
        }
      }

      if (done)
        return *this;

      cur_val = math::zero<val_type>();
      end = (cur_col + 1) * BlockSize;
      for (int i = 0; i < BlockSize; ++i) {
        for (; base[i] && static_cast<ptrdiff_t>(base[i].col()) < end; ++base[i]) {
          cur_val(i, base[i].col() % BlockSize) = base[i].value();
        }
      }

      return *this;
    }

    col_type col() const
    {
      return cur_col;
    }

    val_type value() const
    {
      return cur_val;
    }
  };

  row_iterator row_begin(size_t i) const
  {
    return row_iterator(A, i);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Convert scalar-valued matrix to a block-valued one.
template <class BlockType, class Matrix>
block_matrix_adapter<Matrix, BlockType> block_matrix(const Matrix& A)
{
  return block_matrix_adapter<Matrix, BlockType>(A);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Matrix>
std::shared_ptr<backend::crs<typename math::element_of<
                             typename backend::value_type<Matrix>::type>::type,
                             typename backend::col_type<Matrix>::type,
                             typename backend::ptr_type<Matrix>::type>>
unblock_matrix(const Matrix& B)
{
  typedef typename backend::value_type<Matrix>::type Block;
  typedef typename math::element_of<Block>::type Scalar;
  typedef typename backend::col_type<Matrix>::type Col;
  typedef typename backend::ptr_type<Matrix>::type Ptr;

  const int brows = math::static_rows<Block>::value;
  const int bcols = math::static_cols<Block>::value;

  static_assert(brows > 1 || bcols > 1, "Can not unblock scalar matrix!");

  auto A = std::make_shared<backend::crs<Scalar, Col, Ptr>>();

  A->set_size(backend::rows(B) * brows, backend::cols(B) * bcols);
  A->ptr[0] = 0;

  const ptrdiff_t nb = backend::rows(B);

#pragma omp for
  for (ptrdiff_t ib = 0; ib < nb; ++ib) {
    auto w = backend::row_nonzeros(B, ib);
    for (ptrdiff_t i = 0, ia = ib * brows; i < brows; ++i, ++ia) {
      A->ptr[ia + 1] = w * bcols;
    }
  }

  A->scan_row_sizes();
  A->set_nonzeros();

#pragma omp for
  for (ptrdiff_t ib = 0; ib < nb; ++ib) {
    for (auto b = backend::row_begin(B, ib); b; ++b) {
      auto c = b.col();
      auto v = b.value();

      for (ptrdiff_t i = 0, ia = ib * brows; i < brows; ++i, ++ia) {
        auto row_head = A->ptr[ia];
        for (int j = 0; j < bcols; ++j) {
          A->col[row_head] = c * bcols + j;
          A->val[row_head] = v(i, j);
          ++row_head;
        }
        A->ptr[ia] = row_head;
      }
    }
  }

  std::rotate(A->ptr, A->ptr + A->nrows, A->ptr + A->nrows + 1);
  A->ptr[0] = 0;

  return A;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Matrix>
struct complex_adapter
{
  static_assert(is_complex<typename backend::value_type<Matrix>::type>::value,
                "value type should be complex");

  typedef typename backend::value_type<Matrix>::type::value_type value_type;

  const Matrix& A;

  complex_adapter(const Matrix& A)
  : A(A)
  {}

  size_t rows() const
  {
    return 2 * backend::rows(A);
  }

  size_t cols() const
  {
    return 2 * backend::cols(A);
  }

  size_t nonzeros() const
  {
    return 4 * backend::nonzeros(A);
  }

  struct row_iterator
  {
    typedef typename backend::row_iterator<Matrix>::type Base;
    typedef typename Base::col_type col_type;

    row_iterator(const Base& base, bool row_real)
    : base(base)
    , row_real(row_real)
    , col_real(true)
    {}

    operator bool() const
    {
      return static_cast<bool>(base);
    }

    row_iterator& operator++()
    {
      col_real = !col_real;
      if (col_real)
        ++base;

      return *this;
    }

    col_type col() const
    {
      if (col_real)
        return base.col() * 2;
      else
        return base.col() * 2 + 1;
    }

    value_type value() const
    {
      if (row_real) {
        if (col_real)
          return std::real(base.value());
        else
          return -std::imag(base.value());
      }
      else {
        if (col_real)
          return std::imag(base.value());
        else
          return std::real(base.value());
      }
    }

   private:

    Base base;
    bool row_real;
    bool col_real;
  };

  row_iterator row_begin(size_t i) const
  {
    return row_iterator(backend::row_begin(A, i / 2), i % 2 == 0);
  }
};

template <class Matrix>
complex_adapter<Matrix> complex_matrix(const Matrix& A)
{
  return complex_adapter<Matrix>(A);
}

template <class Range>
boost::iterator_range<typename std::add_pointer<
typename std::conditional<std::is_const<Range>::value,
                          typename std::add_const<
                          typename boost::range_value<
                          typename std::decay<Range>::type>::type::value_type>::type,
                          typename boost::range_value<
                          typename std::decay<Range>::type>::type::value_type>::type>::type>
complex_range(Range& rng)
{
  typedef
  typename std::add_pointer<
  typename std::conditional<
  std::is_const<Range>::value,
  typename std::add_const<
  typename boost::range_value<
  typename std::decay<Range>::type>::type::value_type>::type,
  typename boost::range_value<
  typename std::decay<Range>::type>::type::value_type>::type>::type
  pointer_type;

  pointer_type b = reinterpret_cast<pointer_type>(&rng[0]);
  pointer_type e = b + 2 * boost::size(rng);

  return boost::iterator_range<pointer_type>(b, e);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Generates matrix rows as needed with help of user-provided functor.
/**
 * The generated rows are not stored anywhere.
 */
template <class RowBuilder>
struct matrix_builder
{
  typedef typename RowBuilder::val_type value_type;
  typedef typename RowBuilder::col_type col_type;

  RowBuilder build_row;

  matrix_builder(const RowBuilder& row_builder)
  : build_row(row_builder)
  {}

  size_t rows() const { return build_row.rows(); }
  size_t cols() const { return build_row.rows(); }
  size_t nonzeros() const { return build_row.nonzeros(); }

  struct row_iterator
  {
    typedef typename RowBuilder::col_type col_type;
    typedef typename RowBuilder::val_type val_type;

    typedef typename std::vector<col_type>::const_iterator col_iterator;
    typedef typename std::vector<val_type>::const_iterator val_iterator;

    row_iterator(const RowBuilder& build_row, size_t i)
    : ptr(0)
    {
      build_row(i, m_col, m_val);
    }

    operator bool() const
    {
      return m_col.size() - ptr;
    }

    row_iterator& operator++()
    {
      ++ptr;
      return *this;
    }

    col_type col() const
    {
      return m_col[ptr];
    }

    val_type value() const
    {
      return m_val[ptr];
    }

   private:

    int ptr;
    std::vector<col_type> m_col;
    std::vector<value_type> m_val;
  };

  row_iterator row_begin(size_t i) const
  {
    return row_iterator(build_row, i);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Convenience function returning an instance of matrix_builder<RowBuilder>
template <class RowBuilder>
matrix_builder<RowBuilder> make_matrix(const RowBuilder& row_builder)
{
  return matrix_builder<RowBuilder>(row_builder);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Matrix>
struct reordered_matrix
{
  typedef typename backend::value_type<Matrix>::type value_type;
  typedef typename backend::row_iterator<Matrix>::type base_iterator;

  const Matrix& A;
  const ptrdiff_t* perm;
  const ptrdiff_t* iperm;

  reordered_matrix(const Matrix& A, const ptrdiff_t* perm, const ptrdiff_t* iperm)
  : A(A)
  , perm(perm)
  , iperm(iperm)
  {}

  size_t rows() const
  {
    return backend::rows(A);
  }

  size_t cols() const
  {
    return backend::cols(A);
  }

  size_t nonzeros() const
  {
    return backend::nonzeros(A);
  }

  struct row_iterator
  {
    base_iterator base;
    const ptrdiff_t* iperm;

    row_iterator(const base_iterator& base, const ptrdiff_t* iperm)
    : base(base)
    , iperm(iperm)
    {}

    operator bool() const
    {
      return base;
    }

    row_iterator& operator++()
    {
      ++base;
      return *this;
    }

    ptrdiff_t col() const
    {
      return iperm[base.col()];
    }

    value_type value() const
    {
      return base.value();
    }
  };

  row_iterator row_begin(size_t i) const
  {
    return row_iterator(backend::row_begin(A, perm[i]), iperm);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Vector>
struct reordered_vector
{
  typedef typename backend::value_type<typename std::decay<Vector>::type>::type raw_value_type;
  typedef typename std::conditional<std::is_const<Vector>::value,
                                    const raw_value_type,
                                    raw_value_type>::type value_type;

  Vector& x;
  const ptrdiff_t* perm;

  reordered_vector(Vector& x, const ptrdiff_t* perm)
  : x(x)
  , perm(perm)
  {}

  size_t size() const
  {
    return boost::size(x);
  }

  value_type& operator[](size_t i) const
  {
    return x[perm[i]];
  }

  boost::permutation_iterator<typename std::decay<Vector>::type::iterator,
                              const ptrdiff_t*>
  begin()
  {
    return boost::make_permutation_iterator(boost::begin(x), perm);
  }

  boost::permutation_iterator<typename std::decay<Vector>::type::const_iterator,
                              const ptrdiff_t*>
  begin() const
  {
    return boost::make_permutation_iterator(boost::begin(x), perm);
  }

  boost::permutation_iterator<typename std::decay<Vector>::type::iterator,
                              const ptrdiff_t*>
  end()
  {
    return boost::make_permutation_iterator(boost::end(x), perm + size());
  }

  boost::permutation_iterator<typename std::decay<Vector>::type::const_iterator,
                              const ptrdiff_t*>
  end() const
  {
    return boost::make_permutation_iterator(boost::end(x), perm + size());
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class ordering = Alina::reorder::cuthill_mckee<false>>
class reorder
{
 public:

  template <class Matrix>
  reorder(const Matrix& A)
  : n(backend::rows(A))
  , perm(n)
  , iperm(n)
  {
    ordering::get(A, perm);
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n; ++i)
      iperm[perm[i]] = i;
  }

  template <class Matrix>
  typename std::enable_if<!backend::is_builtin_vector<Matrix>::value,
                          reordered_matrix<Matrix>>::type
  operator()(const Matrix& A) const
  {
    return reordered_matrix<Matrix>(A, perm.data(), iperm.data());
  }

  template <class Vector>
  typename std::enable_if<backend::is_builtin_vector<Vector>::value,
                          reordered_vector<Vector>>::type
  operator()(Vector& x) const
  {
    return reordered_vector<Vector>(x, perm.data());
  }

  template <class Vector>
  typename std::enable_if<backend::is_builtin_vector<Vector>::value,
                          reordered_vector<const Vector>>::type
  operator()(const Vector& x) const
  {
    return reordered_vector<const Vector>(x, perm.data());
  }

  template <class Vector1, class Vector2>
  void forward(const Vector1& x, Vector2& y) const
  {
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n; ++i)
      y[i] = x[perm[i]];
  }

  template <class Vector1, class Vector2>
  void inverse(const Vector1& x, Vector2& y) const
  {
#pragma omp parallel for
    for (ptrdiff_t i = 0; i < n; ++i)
      y[perm[i]] = x[i];
  }

 private:

  ptrdiff_t n;
  backend::numa_vector<ptrdiff_t> perm, iperm;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Matrix, class Scale>
struct scaled_matrix
{
  typedef typename backend::value_type<Matrix>::type value_type;
  typedef typename backend::value_type<Scale>::type scale_type;

  const Matrix& A;
  const Scale& s;

  scaled_matrix(const Matrix& A, const Scale& s)
  : A(A)
  , s(s)
  {}

  size_t rows() const { return backend::rows(A); }
  size_t cols() const { return backend::cols(A); }
  size_t nonzeros() const { return backend::nonzeros(A); }

  struct row_iterator : public backend::row_iterator<Matrix>::type
  {
    typedef typename backend::row_iterator<Matrix>::type Base;

    scale_type si;
    const Scale& s;

    row_iterator(const Matrix& A, const Scale& s, size_t i)
    : Base(A, i)
    , si(s[i])
    , s(s)
    {}

    value_type value() const
    {
      return si * static_cast<const Base*>(this)->value() * s[this->col()];
    }
  };

  row_iterator row_begin(size_t i) const
  {
    return row_iterator(A, s, i);
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class Scale>
struct scaled_problem
{
  typedef typename Backend::params backend_params;

  const std::shared_ptr<Scale> s;
  const backend_params& bprm;

  scaled_problem(std::shared_ptr<Scale> s, const backend_params& bprm = backend_params())
  : s(s)
  , bprm(bprm)
  {}

  template <class Matrix>
  scaled_matrix<Matrix, Scale> matrix(const Matrix& A) const
  {
    return scaled_matrix<Matrix, Scale>(A, *s);
  }

  template <class Vector>
  std::shared_ptr<typename Backend::vector> rhs(const Vector& v) const
  {
    auto t = Backend::copy_vector(v, bprm);
    (*this)(*t);
    return t;
  }

  template <class Vector>
  void operator()(Vector& x) const
  {
    typedef typename backend::value_type<Vector>::type value_type;
    typedef typename math::scalar_of<value_type>::type scalar_type;

    const auto one = math::identity<scalar_type>();
    const auto zero = math::zero<scalar_type>();

    if (backend::is_builtin_vector<Vector>::value) {
      backend::vmul(one, *s, x, zero, x);
    }
    else {
      backend::vmul(one, *Backend::copy_vector(*s, bprm), x, zero, x);
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class Matrix>
scaled_problem<Backend,
               std::vector<
               typename math::scalar_of<
               typename backend::value_type<Matrix>::type>::type>>
scale_diagonal(
const Matrix& A,
const typename Backend::params& bprm = typename Backend::params())
{
  typedef typename backend::value_type<Matrix>::type value_type;
  typedef typename math::scalar_of<value_type>::type scalar_type;
  ptrdiff_t n = backend::rows(A);
  auto s = std::make_shared<std::vector<scalar_type>>(n);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < n; ++i) {
    for (auto a = backend::row_begin(A, i); a; ++a) {
      if (a.col() == i) {
        (*s)[i] = math::inverse(sqrt(math::norm(a.value())));
        break;
      }
    }
  }

  return scaled_problem<Backend, std::vector<scalar_type>>(s, bprm);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <typename Ptr, typename Col, typename Val>
std::shared_ptr<backend::crs<Val>>
zero_copy(size_t nrows, size_t ncols, const Ptr* ptr, const Col* col, const Val* val)
{
  // Check that Ptr and Col types are binary-compatible with ptrdiff_t:
  static_assert(std::is_integral<Ptr>::value, "Unsupported Ptr type");
  static_assert(std::is_integral<Col>::value, "Unsupported Col type");
  static_assert(sizeof(Ptr) == sizeof(ptrdiff_t), "Unsupported Ptr type");
  static_assert(sizeof(Col) == sizeof(ptrdiff_t), "Unsupported Col type");

  auto A = std::make_shared<backend::crs<Val>>();
  A->nrows = nrows;
  A->ncols = ncols;
  A->nnz = nrows ? ptr[nrows] : 0;

  A->ptr = (ptrdiff_t*)ptr;
  A->col = (ptrdiff_t*)col;
  A->val = (Val*)val;

  A->own_data = false;

  return A;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <typename Ptr, typename Col, typename Val>
std::shared_ptr<backend::crs<Val>>
zero_copy(size_t n, const Ptr* ptr, const Col* col, const Val* val)
{
  return zero_copy(n, n, ptr, col, val);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <typename Ptr, typename Col, typename Val>
std::shared_ptr<backend::crs<Val, Col, Ptr>>
zero_copy_direct(size_t nrows, size_t ncols, const Ptr* ptr, const Col* col, const Val* val)
{
  auto A = std::make_shared<backend::crs<Val, Col, Ptr>>();
  A->nrows = nrows;
  A->ncols = ncols;
  A->nnz = nrows ? ptr[nrows] : 0;

  A->ptr = const_cast<Ptr*>(ptr);
  A->col = const_cast<Col*>(col);
  A->val = const_cast<Val*>(val);

  A->own_data = false;

  return A;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <typename Ptr, typename Col, typename Val>
std::shared_ptr<backend::crs<Val, Col, Ptr>>
zero_copy_direct(size_t n, const Ptr* ptr, const Col* col, const Val* val)
{
  return zero_copy_direct(n, n, ptr, col, val);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::adapter

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::backend
{
template <class Vector>
struct is_builtin_vector<adapter::reordered_vector<Vector>>
: is_builtin_vector<typename std::decay<Vector>::type>
{};
} // namespace Arcane::Alina::backend

namespace Arcane::Alina::backend::detail
{

template <class Matrix, class BlockType>
struct use_builtin_matrix_ops<adapter::block_matrix_adapter<Matrix, BlockType>>
: std::true_type
{};

template <class Matrix>
struct use_builtin_matrix_ops<Alina::adapter::complex_adapter<Matrix>>
: std::true_type
{};

template <class RowBuilder>
struct use_builtin_matrix_ops<Alina::adapter::matrix_builder<RowBuilder>>
: std::true_type
{};

template <typename N, typename PRng, typename CRng, typename VRng>
struct use_builtin_matrix_ops<std::tuple<N, PRng, CRng, VRng>>
: std::true_type
{};

template <class Matrix>
struct use_builtin_matrix_ops<adapter::reordered_matrix<Matrix>>
: std::true_type
{};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::backend::detail

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
