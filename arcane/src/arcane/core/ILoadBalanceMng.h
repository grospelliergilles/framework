﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ILoadBalanceMng.h                                           (C) 2000-2023 */
/*                                                                           */
/* Interface de description des caracteristiques du probleme pour le module  */
/* d'equilibrage de charge.                                                  */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ILOADBALANCEMNG_H
#define ARCANE_ILOADBALANCEMNG_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/ArcaneTypes.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Interface d'enregistrement des variables pour l'equilibrage de charge.
 */
class ILoadBalanceMng
{
 public:

  virtual ~ILoadBalanceMng() {} //!< Libère les ressources.

 public:

  /*!
   * Methodes utilisees par les modules clients pour definir les criteres
   * de partitionnement.
   */
  virtual void addMass(VariableCellInt32& count, const String& entity="") =0;
  virtual void addCriterion(VariableCellInt32& count) =0;
  virtual void addCriterion(VariableCellReal& count) =0;
  virtual void addCommCost(VariableFaceInt32& count, const String& entity="") =0;

  virtual void reset() =0;

  /*!
   * Methodes utilisees par le MeshPartitioner pour acceder a la description
   * du probleme.
   */
  virtual void setMassAsCriterion(bool active=true) =0;
  virtual void setNbCellsAsCriterion(bool active=true) =0;
  virtual Integer nbCriteria() =0;
  virtual void setCellCommContrib(bool active=true) =0;
  virtual bool cellCommContrib() const =0;
  virtual void setComputeComm(bool active=true) =0;
  virtual void initAccess(IMesh* mesh) =0;
  virtual const VariableFaceReal& commCost() const =0;
  virtual const VariableCellReal& massWeight() const =0;
  virtual const VariableCellReal& massResWeight() const =0;
  virtual const VariableCellArrayReal& mCriteriaWeight() const =0;
  virtual void endAccess() =0;
  virtual void notifyEndPartition() =0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // End namespace Arcane

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
