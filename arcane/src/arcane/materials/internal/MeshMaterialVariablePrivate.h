﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* MeshMaterialVariablePrivate.h                               (C) 2000-2023 */
/*                                                                           */
/* Partie privée d'une variable sur un matériau du maillage.                 */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_MATERIALS_INTERNAL_MESHMATERIALVARIABLEPRIVATE_H
#define ARCANE_MATERIALS_INTERNAL_MESHMATERIALVARIABLEPRIVATE_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/Array.h"
#include "arcane/utils/ScopedPtr.h"

#include "arcane/core/VariableDependInfo.h"
#include "arcane/core/materials/internal/IMeshMaterialVariableInternal.h"

#include "arcane/materials/MeshMaterialVariableDependInfo.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane
{
class IObserver;
}

namespace Arcane::Materials
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \internal
 * \brief Partie privée d'une variable matériau.
 */
class MeshMaterialVariablePrivate
: public IMeshMaterialVariableInternal
{
 public:

  MeshMaterialVariablePrivate(const MaterialVariableBuildInfo& v, MatVarSpace mvs,
                              MeshMaterialVariable* variable);
  ~MeshMaterialVariablePrivate();

 public:

  MatVarSpace space() const { return m_var_space; }
  bool hasRecursiveDepend() const { return m_has_recursive_depend; }
  const String& name() const { return m_name; }
  IMeshMaterialMng* materialMng() const { return m_material_mng; }
  IMeshMaterialVariableInternal* _internalApi() { return this; }

 public:

  Int32 dataTypeSize() const override;
  void copyToBuffer(SmallSpan<const MatVarIndex> matvar_indexes,
                    Span<std::byte> bytes, RunQueue* queue) const override;

  void copyFromBuffer(SmallSpan<const MatVarIndex> matvar_indexes,
                      Span<const std::byte> bytes, RunQueue* queue) override;

  Ref<IData> internalCreateSaveDataRef(Integer nb_value) override;

  void saveData(IMeshComponent* component,IData* data) override;

  void restoreData(IMeshComponent* component,IData* data,Integer data_index,Int32ConstArrayView ids,bool allow_null_id) override;

  void copyGlobalToPartial(Int32 var_index,Int32ConstArrayView local_ids,Int32ConstArrayView indexes_in_multiple) override;

  void copyPartialToGlobal(Int32 var_index,Int32ConstArrayView local_ids,Int32ConstArrayView indexes_in_multiple) override;

  void initializeNewItems(const ComponentItemListBuilder& list_builder) override;

 public:

  Int32 m_nb_reference;
  MeshMaterialVariableRef* m_first_reference; //! Première référence sur la variable

 private:

  String m_name;
  IMeshMaterialMng* m_material_mng;

 public:

  /*!
   * \brief Stocke les références sur les variables tableaux qui servent pour
   * stocker les valeurs par matériau.
   * Il faut garder une référence pour éviter que la variable ne soit détruite
   * si elle n'est plus utilisée par ailleurs.
   */
  UniqueArray<VariableRef*> m_refs;

  bool m_keep_on_change;
  IObserver* m_global_variable_changed_observer;

  //! Liste des dépendances de cette variable
  UniqueArray<MeshMaterialVariableDependInfo> m_mat_depends;

  //! Liste des dépendances de cette variable
  UniqueArray<VariableDependInfo> m_depends;

  //! Tag de la dernière modification par matériau
  UniqueArray<Int64> m_modified_times;

  //! Fonction de calcul
  ScopedPtrT<IMeshMaterialVariableComputeFunction> m_compute_function;

 private:

  bool m_has_recursive_depend;
  MatVarSpace m_var_space;
  MeshMaterialVariable* m_variable = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Materials

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

