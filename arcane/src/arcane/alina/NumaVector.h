// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* NumaVector.h                                                (C) 2026-2026 */
/*                                                                           */
/* NUMA-aware vector container.                                              */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_NUMAVECTOR_H
#define ARCANE_ALINA_NUMAVECTOR_H
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "arcane/alina/AlinaUtils.h"
#include "arcane/alina/ValueTypeInterface.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::backend
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief NUMA-aware vector container.
 */
template <class T>
class numa_vector
{
 public:

  typedef T value_type;

  numa_vector()
  : n(0)
  , p(0)
  {}

  explicit numa_vector(size_t n, bool init = true)
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
    static_assert(std::is_same<std::random_access_iterator_tag,
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

} // namespace Arcane::Alina::backend

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
