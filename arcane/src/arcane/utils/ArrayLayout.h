﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArrayLayout.h                                               (C) 2000-2022 */
/*                                                                           */
/* Gestion de la disposition mémoire pour les tableaux N-dimensions.         */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_UTILS_ARRAYLAYOUT_H
#define ARCANE_UTILS_ARRAYLAYOUT_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/ArrayBoundsIndex.h"

#include <array>

/*
 * ATTENTION:
 *
 * Toutes les classes de ce fichier sont expérimentales et l'API n'est pas
 * figée. A NE PAS UTILISER EN DEHORS DE ARCANE.
 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<Int32 I,Int32 J>
class ArrayLayout2
{
 public:

  static constexpr Int64 LastExtent = J;

  static ARCCORE_HOST_DEVICE constexpr Int64
  offset(ArrayBoundsIndex<2> idx,Int64 extent1)
  {
    return (extent1 * idx[I]) + idx[J];
  }

  static std::array<Int32,2> layoutInfo() { return { I, J }; }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<Int32 I,Int32 J,Int32 K>
class ArrayLayout3
{
 public:

  static constexpr Int64 LastExtent = K;

  static ARCCORE_HOST_DEVICE constexpr Int64
  offset(ArrayBoundsIndex<3> idx,Int64 extent1,Int64 extent2)
  {
    return (extent2 * idx[I]) + (extent1*idx[J]) + idx[K];
  }

  static ARCCORE_HOST_DEVICE constexpr Int64
  computeOffsetIndexes(std::array<Int64,3> extents)
  {
    return extents[J] * extents[K];
  }

  static std::array<Int32,3> layoutInfo() { return { I, J, K }; }
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Layout par défaut pour chaque dimension

template<> class DefaultLayout<2> : public ArrayLayout2<0,1> {};
template<> class DefaultLayout<3> : public ArrayLayout3<0,1,2> {};
//template<> class DefaultLayout<3> : public ArrayLayout3<2,1,0> {};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // End namespace Arcane

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
