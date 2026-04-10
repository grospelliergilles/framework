// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedCoarseningRuntime.h                              (C) 2026-2026 */
/*                                                                           */
/* Runtime wrapper for distributed coarsening schemes.                       */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MPI_DISTRIBUTEDCOARSENINGRUNTIME_H
#define ARCANE_ALINA_MPI_DISTRIBUTEDCOARSENINGRUNTIME_H
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

#include <arcane/alina/mpi/DistributedMatrix.h>
#include <arcane/alina/mpi/DistributedCoarsening.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::mpi::coarsening
{

enum type
{
  aggregation,
  smoothed_aggregation
};

inline std::ostream& operator<<(std::ostream& os, type s)
{
  switch (s) {
  case aggregation:
    return os << "aggregation";
  case smoothed_aggregation:
    return os << "smoothed_aggregation";
  default:
    return os << "???";
  }
}

inline std::istream& operator>>(std::istream& in, type& s)
{
  std::string val;
  in >> val;

  if (val == "aggregation")
    s = aggregation;
  else if (val == "smoothed_aggregation")
    s = smoothed_aggregation;
  else
    throw std::invalid_argument("Invalid coarsening value. Valid choices are: "
                                "aggregation, smoothed_aggregation.");

  return in;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct DistributedCoarseningRuntime
{
  typedef DistributedMatrix<Backend> matrix;
  typedef PropertyTree params;

  type c;
  void* handle = nullptr;

  explicit DistributedCoarseningRuntime(params prm = params())
  : c(prm.get("type", smoothed_aggregation))
  {
    if (!prm.erase("type"))
      ARCANE_ALINA_PARAM_MISSING("type");

    switch (c) {
    case aggregation: {
      typedef Alina::mpi::coarsening::aggregation<Backend> C;
      handle = static_cast<void*>(new C(prm));
    } break;
    case smoothed_aggregation: {
      typedef Alina::mpi::coarsening::smoothed_aggregation<Backend> C;
      handle = static_cast<void*>(new C(prm));
    } break;
    default:
      throw std::invalid_argument("Unsupported coarsening type");
    }
  }

  ~DistributedCoarseningRuntime()
  {
    switch (c) {
    case aggregation: {
      typedef Alina::mpi::coarsening::aggregation<Backend> C;
      delete static_cast<C*>(handle);
    } break;
    case smoothed_aggregation: {
      typedef Alina::mpi::coarsening::smoothed_aggregation<Backend> C;
      delete static_cast<C*>(handle);
    } break;
    default:
      break;
    }
  }

  std::tuple<std::shared_ptr<matrix>, std::shared_ptr<matrix>>
  transfer_operators(const matrix& A)
  {
    switch (c) {
    case aggregation: {
      typedef Alina::mpi::coarsening::aggregation<Backend> C;
      return static_cast<C*>(handle)->transfer_operators(A);
    }
    case smoothed_aggregation: {
      typedef Alina::mpi::coarsening::smoothed_aggregation<Backend> C;
      return static_cast<C*>(handle)->transfer_operators(A);
    }
    default:
      throw std::invalid_argument("Unsupported partition type");
    }
  }

  std::shared_ptr<matrix>
  coarse_operator(const matrix& A, const matrix& P, const matrix& R) const
  {
    switch (c) {
    case aggregation: {
      typedef Alina::mpi::coarsening::aggregation<Backend> C;
      return static_cast<C*>(handle)->coarse_operator(A, P, R);
    }
    case smoothed_aggregation: {
      typedef Alina::mpi::coarsening::smoothed_aggregation<Backend> C;
      return static_cast<C*>(handle)->coarse_operator(A, P, R);
    }
    default:
      throw std::invalid_argument("Unsupported partition type");
    }
  }
};

template <class Backend>
unsigned block_size(const DistributedCoarseningRuntime<Backend>& w)
{
  switch (w.c) {
  case aggregation: {
    typedef Alina::mpi::coarsening::aggregation<Backend> C;
    return block_size(*static_cast<const C*>(w.handle));
  }
  case smoothed_aggregation: {
    typedef Alina::mpi::coarsening::smoothed_aggregation<Backend> C;
    return block_size(*static_cast<const C*>(w.handle));
  }
  default:
    throw std::invalid_argument("Unsupported coarsening type");
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::runtime::mpi::coarsening

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
