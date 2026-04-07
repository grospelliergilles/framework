// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* mp_partition_runtime.h                                      (C) 2026-2026 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MPI_MP_PARTITION_RUNTIME_H
#define ARCANE_ALINA_MPI_MP_PARTITION_RUNTIME_H
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

#include <memory>

#include <arcane/alina/util.h>
#include <arcane/alina/mpi/mp_partition_merge.h>
#if defined(ARCANE_ALINA_HAVE_PARMETIS)
#include <arcane/alina/mpi/mp_partition_parmetis.h>
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::mpi::partition
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

enum type
{
  merge
#ifdef ARCANE_ALINA_HAVE_PARMETIS
  ,
  parmetis
#endif
};

inline std::ostream&
operator<<(std::ostream& os, type s)
{
  switch (s) {
  case merge:
    return os << "merge";
#ifdef ARCANE_ALINA_HAVE_PARMETIS
  case parmetis:
    return os << "parmetis";
#endif
  default:
    return os << "???";
  }
}

inline std::istream&
operator>>(std::istream& in, type& s)
{
  std::string val;
  in >> val;

  if (val == "merge")
    s = merge;
#ifdef ARCANE_ALINA_HAVE_PARMETIS
  else if (val == "parmetis")
    s = parmetis;
#endif
  else
    throw std::invalid_argument("Invalid partitioner value. Valid choices are: "
                                "merge"
#ifdef ARCANE_ALINA_HAVE_PARMETIS
                                ", parmetis"
#endif
                                ".");

  return in;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend>
struct wrapper
{
  typedef Alina::mpi::distributed_matrix<Backend> matrix;
  typedef Alina::PropertyTree params;

  type t;
  void* handle;

  wrapper(params prm = params())
  : t(prm.get("type",
#if defined(ARCANE_ALINA_HAVE_PARMETIS)
              parmetis
#else
              merge
#endif
              ))
  , handle(0)
  {
    if (!prm.erase("type"))
      ARCANE_ALINA_PARAM_MISSING("type");

    switch (t) {
    case merge: {
      typedef Alina::mpi::partition::merge<Backend> R;
      handle = static_cast<void*>(new R(prm));
    } break;
#ifdef ARCANE_ALINA_HAVE_PARMETIS
    case parmetis: {
      typedef Alina::mpi::partition::parmetis<Backend> R;
      handle = static_cast<void*>(new R(prm));
    } break;
#endif
    default:
      throw std::invalid_argument("Unsupported partition type");
    }
  }

  ~wrapper()
  {
    switch (t) {
    case merge: {
      typedef Alina::mpi::partition::merge<Backend> R;
      delete static_cast<R*>(handle);
    } break;
#ifdef ARCANE_ALINA_HAVE_PARMETIS
    case parmetis: {
      typedef Alina::mpi::partition::parmetis<Backend> R;
      delete static_cast<R*>(handle);
    } break;
#endif
    default:
      break;
    }
  }

  bool is_needed(const matrix& A) const
  {
    switch (t) {
    case merge: {
      typedef Alina::mpi::partition::merge<Backend> R;
      return static_cast<const R*>(handle)->is_needed(A);
    }
#ifdef ARCANE_ALINA_HAVE_PARMETIS
    case parmetis: {
      typedef Alina::mpi::partition::parmetis<Backend> R;
      return static_cast<const R*>(handle)->is_needed(A);
    }
#endif
    default:
      throw std::invalid_argument("Unsupported partition type");
    }
  }

  std::shared_ptr<matrix> operator()(const matrix& A, unsigned block_size = 1) const
  {
    switch (t) {
    case merge: {
      typedef Alina::mpi::partition::merge<Backend> R;
      return static_cast<const R*>(handle)->operator()(A, block_size);
    }
#ifdef ARCANE_ALINA_HAVE_PARMETIS
    case parmetis: {
      typedef Alina::mpi::partition::parmetis<Backend> R;
      return static_cast<const R*>(handle)->operator()(A, block_size);
    }
#endif
    default:
      throw std::invalid_argument("Unsupported partition type");
    }
  }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::runtime::mpi::partition

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
