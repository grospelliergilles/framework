// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* DistributedSolver.h                                         (C) 2026-2026 */
/*                                                                           */
/* Adapters to handle distribution other standard solvers.                   */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_MPI_DISTRIBUTEDSOLVER_H
#define ARCANE_ALINA_MPI_DISTRIBUTEDSOLVER_H
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

#include <arcane/alina/ConjugateGradientSolver.h>
#include <arcane/alina/BiCGStabSolver.h>
#include <arcane/alina/BiCGStabLSolver.h>
#include <arcane/alina/FlexibleGMRESSolver.h>
#include <arcane/alina/GMRESSolver.h>
#include <arcane/alina/IDRSSolver.h>
#include <arcane/alina/LooseGMRESSolver.h>
#include <arcane/alina/PreconditionerOnlySolver.h>
#include <arcane/alina/RichardsonSolver.h>
#include <arcane/alina/DistributedInnerProduct.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedConjugateGradientSolver
: public Alina::solver::ConjugateGradientSolver<Backend, InnerProduct>
{
  typedef Alina::solver::ConjugateGradientSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedBiCGStabSolver
: public Alina::solver::BiCGStabSolver<Backend, InnerProduct>
{
  typedef Alina::solver::BiCGStabSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedBiCGStabLSolver
: public Alina::solver::BiCGStabLSolver<Backend, InnerProduct>
{
  typedef Alina::solver::BiCGStabLSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedFlexibleGMRESSolver
: public Alina::solver::FlexibleGMRESSolver<Backend, InnerProduct>
{
  typedef Alina::solver::FlexibleGMRESSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedGMRESSolver
: public Alina::solver::GMRESSolver<Backend, InnerProduct>
{
  typedef Alina::solver::GMRESSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedIDRSSolver
: public Alina::solver::IDRSSolver<Backend, InnerProduct>
{
  typedef Alina::solver::IDRSSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedLooseGMRESSolver
: public Alina::solver::LooseGMRESSolver<Backend, InnerProduct>
{
  typedef Alina::solver::LooseGMRESSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class DistributedPreconditionerOnlySolver
: public Alina::solver::PreconditionerOnlySolver<Backend, InnerProduct>
{
  typedef Alina::solver::PreconditionerOnlySolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <class Backend, class InnerProduct = DistributedInnerProduct>
class RichardsonSolve
: public Alina::solver::RichardsonSolver<Backend, InnerProduct>
{
  typedef Alina::solver::RichardsonSolver<Backend, InnerProduct> Base;

 public:

  using Base::Base;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
