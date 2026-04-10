// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedSolverRuntime.h                                  (C) 2000-2026 */
/*                                                                           */
/* Runtime-configurable MPI wrapper around iterative solvers.                */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_DISTRIBUTEDSOLVERRUNTIME_H
#define ARCANE_ALINA_DISTRIBUTEDSOLVERRUNTIME_H
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

#include <arcane/alina/SolverRuntime.h>
#include <arcane/alina/mpi/DistributedInnerProduct.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime::mpi::solver
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = Alina::mpi::inner_product>
struct DistributedSolverRuntime
: public Alina::runtime::solver::SolverRuntime<Backend, InnerProduct>
{
  typedef Alina::runtime::solver::SolverRuntime<Backend, InnerProduct> Base;
  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::runtime::mpi::solver

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
