// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* mp_solver.h                                                 (C) 2026-2026 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MPI_MP_SOLVER_H
#define ARCANE_ALINA_MPI_MP_SOLVER_H
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

#include <arcane/alina/solver_cg.h>
#include <arcane/alina/solver_bicgstab.h>
#include <arcane/alina/solver_bicgstabl.h>
#include <arcane/alina/solver_fgmres.h>
#include <arcane/alina/solver_gmres.h>
#include <arcane/alina/solver_idrs.h>
#include <arcane/alina/solver_lgmres.h>
#include <arcane/alina/solver_preonly.h>
#include <arcane/alina/solver_richardson.h>
#include <arcane/alina/mpi/mp_inner_product.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::mpi::solver
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class cg
: public Alina::solver::ConjugateGradient<Backend, InnerProduct>
{
  typedef Alina::solver::ConjugateGradient<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class bicgstab
: public Alina::solver::BiCGStabSolver<Backend, InnerProduct>
{
  typedef Alina::solver::BiCGStabSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class bicgstabl
: public Alina::solver::bicgstabl<Backend, InnerProduct>
{
  typedef Alina::solver::bicgstabl<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class fgmres
: public Alina::solver::fgmres<Backend, InnerProduct>
{
  typedef Alina::solver::fgmres<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class gmres
: public Alina::solver::GMRESSolver<Backend, InnerProduct>
{
  typedef Alina::solver::GMRESSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class idrs
: public Alina::solver::idrs<Backend, InnerProduct>
{
  typedef Alina::solver::idrs<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class lgmres
: public Alina::solver::lgmres<Backend, InnerProduct>
{
  typedef Alina::solver::lgmres<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class preonly
: public Alina::solver::preonly<Backend, InnerProduct>
{
  typedef Alina::solver::preonly<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = mpi::inner_product>
class richardson
: public Alina::solver::RichardsonSolver<Backend, InnerProduct>
{
  typedef Alina::solver::RichardsonSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::mpi::solver

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
