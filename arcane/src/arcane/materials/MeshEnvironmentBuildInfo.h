﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MeshEnvironmentBuildInfo.h                                  (C) 2000-2023 */
/*                                                                           */
/* Informations pour la création d'un milieu.                                */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_MATERIALS_MESHENVIRONMENTBUILDINFO_H
#define ARCANE_MATERIALS_MESHENVIRONMENTBUILDINFO_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/String.h"
#include "arcane/utils/Array.h"
#include "arcane/materials/MaterialsGlobal.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Materials
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup ArcaneMaterials
 * \brief Informations pour la création d'un milieu.
 *
 * Cette instance contient les infos nécessaire à la création d'un milieu.
 * Une fois les informations spécifiées, il faut créer le milieu
 * via IMeshMaterialMng::createEnvironment().
 *
 * Pour l'instant, la seule information pertinante sur un milieu est son
 * nom et la liste des matériaux le composant.
 */
class ARCANE_MATERIALS_EXPORT MeshEnvironmentBuildInfo
{
 public:
  class MatInfo
  {
   public:
    MatInfo(const String& name) : m_name(name){}
   public:
    String m_name;
   public:
    // Le constructeur vide ne doit pas être dispo mais ca plante à
    // la compilation avec VS2010 s'il est absent
    MatInfo() {}
  };
 public:

  MeshEnvironmentBuildInfo(const String& name);
  ~MeshEnvironmentBuildInfo();

 public:

  //! Nom du milieu
  const String& name() const { return m_name; }

  /*!
   * \brief Ajoute le matériau de nom \a name au milieu
   *
   * Le matériau doit déjà avoir été enregistré via
   * IMeshMaterialMng::registerMaterialInfo().
   */
  void addMaterial(const String& name);

 public:

  /*!
   * \internal
   * Liste des matériaux.
   */
  ConstArrayView<MatInfo> materials() const
  {
    return m_materials;
  }

 private:

  String m_name;
  UniqueArray<MatInfo> m_materials;

  void _checkValid(const String& name);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

