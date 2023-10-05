// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2023 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* .cc                                 (C) 2000-2023 */
/*                                                                           */
/* .                  */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "CartesianMeshNumberingMng.h"

#include "arcane/core/IMesh.h"
#include "arcane/core/ItemPrinter.h"
#include "arcane/core/IItemFamily.h"
#include "arcane/core/IParallelMng.h"
#include "arcane/core/VariableTypes.h"
#include "arcane/core/Properties.h"
#include "arcane/core/IMeshModifier.h"
#include "arcane/core/MeshStats.h"
#include "arcane/core/ICartesianMeshGenerationInfo.h"
#include "arcane/core/MeshEvents.h"
#include "arcane/utils/Real3.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/


CartesianMeshNumberingMng::
CartesianMeshNumberingMng(IMesh* mesh)
: TraceAccessor(mesh->traceMng())
, m_mesh(mesh)
, m_pattern(2)
{
  auto* m_generation_info = ICartesianMeshGenerationInfo::getReference(m_mesh,true);

  Int64ConstArrayView global_nb_cells_by_direction = m_generation_info->globalNbCells();
  m_nb_cell_x = global_nb_cells_by_direction[MD_DirX];
  if (m_nb_cell_x <= 0)
    ARCANE_FATAL("Bad value '{0}' for globalNbCells()[MD_DirX] (should be >0)", m_nb_cell_x);

  m_nb_cell_y = global_nb_cells_by_direction[MD_DirY];
  if (m_nb_cell_y <= 0)
    ARCANE_FATAL("Bad value '{0}' for globalNbCells()[MD_DirY] (should be >0)", m_nb_cell_y);

  m_nb_cell_z = global_nb_cells_by_direction[MD_DirZ];
  if (m_nb_cell_z <= 0)
    ARCANE_FATAL("Bad value '{0}' for globalNbCells()[MD_DirZ] (should be >0)", m_nb_cell_z);

  m_first_cell_uid_level.add(0);
  m_first_node_uid_level.add(0);
  m_first_face_uid_level.add(0);
}

Int64 CartesianMeshNumberingMng::
getFirstCellUidLevel(Integer level)
{
  if(level < m_first_cell_uid_level.size()){
   return m_first_cell_uid_level[level];
  }

  const Int32 dimension = m_mesh->dimension();

  Int64 level_i_nb_cell_x = getGlobalNbCellsX(m_nb_face_level.size());
  Int64 level_i_nb_cell_y = getGlobalNbCellsY(m_nb_face_level.size());
  Int64 level_i_nb_cell_z = getGlobalNbCellsZ(m_nb_face_level.size());

  for(Integer i = m_nb_cell_level.size(); i < level+1; ++i){
    if(dimension == 2){
      m_nb_cell_level.add(level_i_nb_cell_x * level_i_nb_cell_y);
    }
    else if(dimension == 3){
      m_nb_cell_level.add(level_i_nb_cell_x * level_i_nb_cell_y * level_i_nb_cell_z);
    }
    m_first_cell_uid_level.add(m_first_cell_uid_level[i] + m_nb_cell_level[i]);

    level_i_nb_cell_x *= m_pattern;
    level_i_nb_cell_y *= m_pattern;
    level_i_nb_cell_z *= m_pattern;
  }

  return m_first_cell_uid_level[level];
}


Int64 CartesianMeshNumberingMng::
getFirstNodeUidLevel(Integer level)
{
  if(level < m_first_node_uid_level.size()){
    return m_first_node_uid_level[level];
  }

  const Int32 dimension = m_mesh->dimension();

  Int64 level_i_nb_cell_x = getGlobalNbCellsX(m_nb_face_level.size());
  Int64 level_i_nb_cell_y = getGlobalNbCellsY(m_nb_face_level.size());
  Int64 level_i_nb_cell_z = getGlobalNbCellsZ(m_nb_face_level.size());

  for(Integer i = m_nb_node_level.size(); i < level+1; ++i){
    if(dimension == 2){
      m_nb_node_level.add((level_i_nb_cell_x + 1) * (level_i_nb_cell_y + 1));
    }
    else if(dimension == 3){
      m_nb_node_level.add((level_i_nb_cell_x + 1) * (level_i_nb_cell_y + 1) * (level_i_nb_cell_z + 1));
    }
    m_first_node_uid_level.add(m_first_node_uid_level[i] + m_nb_node_level[i]);

    level_i_nb_cell_x *= m_pattern;
    level_i_nb_cell_y *= m_pattern;
    level_i_nb_cell_z *= m_pattern;
  }

  return m_first_node_uid_level[level];
}


Int64 CartesianMeshNumberingMng::
getFirstFaceUidLevel(Integer level)
{
  if(level < m_first_face_uid_level.size()){
    return m_first_face_uid_level[level];
  }

  const Int32 dimension = m_mesh->dimension();

  Int64 level_i_nb_cell_x = getGlobalNbCellsX(m_nb_face_level.size());
  Int64 level_i_nb_cell_y = getGlobalNbCellsY(m_nb_face_level.size());
  Int64 level_i_nb_cell_z = getGlobalNbCellsZ(m_nb_face_level.size());

  for(Integer i = m_nb_face_level.size(); i < level+1; ++i){
    if(dimension == 2){
      m_nb_face_level.add((level_i_nb_cell_x * level_i_nb_cell_y) * 2 + level_i_nb_cell_x*2 + level_i_nb_cell_y);
    }
    else if(dimension == 3){
      m_nb_face_level.add((level_i_nb_cell_z + 1) * level_i_nb_cell_x * level_i_nb_cell_y
                        + (level_i_nb_cell_x + 1) * level_i_nb_cell_y * level_i_nb_cell_z
                        + (level_i_nb_cell_y + 1) * level_i_nb_cell_z * level_i_nb_cell_x);
    }
    m_first_face_uid_level.add(m_first_face_uid_level[i] + m_nb_face_level[i]);

    level_i_nb_cell_x *= m_pattern;
    level_i_nb_cell_y *= m_pattern;
    level_i_nb_cell_z *= m_pattern;
  }

  return m_first_face_uid_level[level];
}

Int64 CartesianMeshNumberingMng::
getGlobalNbCellsX(Integer level)
{
  return m_nb_cell_x * std::pow(m_pattern, level);
}

Int64 CartesianMeshNumberingMng::
getGlobalNbCellsY(Integer level)
{
  return m_nb_cell_y * std::pow(m_pattern, level);
}

Int64 CartesianMeshNumberingMng::
getGlobalNbCellsZ(Integer level)
{
  return m_nb_cell_z * std::pow(m_pattern, level);
}

Integer CartesianMeshNumberingMng::
getPattern()
{
  return m_pattern;
}

// Tant que l'on a un unique "pattern" pour x, y, z, pas besoin de trois méthodes.
Int64 CartesianMeshNumberingMng::
getOffsetLevelToLevel(Int64 coord, Integer level_from, Integer level_to)
{
  ARCANE_ASSERT((level_from < level_to), ("Pb level_from level_to"));
  return coord * m_pattern * (level_to - level_from);
}

// TODO : Spécialiser pour 2D ?
Int64 CartesianMeshNumberingMng::
uidToCoordX(Int64 uid, Integer level)
{
  Int64 nb_cell_x = getGlobalNbCellsX(level);
  Int64 nb_cell_y = getGlobalNbCellsY(level);

  Int64 to2d = uid % (nb_cell_x * nb_cell_y);
  return to2d % nb_cell_x;
}

// TODO : Spécialiser pour 2D ?
Int64 CartesianMeshNumberingMng::
uidToCoordY(Int64 uid, Integer level)
{
  Int64 nb_cell_x = getGlobalNbCellsX(level);
  Int64 nb_cell_y = getGlobalNbCellsY(level);

  Int64 to2d = uid % (nb_cell_x * nb_cell_y);
  return to2d / nb_cell_x;
}

Int64 CartesianMeshNumberingMng::
uidToCoordZ(Int64 uid, Integer level)
{
  Int64 nb_cell_x = getGlobalNbCellsX(level);
  Int64 nb_cell_y = getGlobalNbCellsY(level);

  return uid / (nb_cell_x * nb_cell_y);
}



Int64 CartesianMeshNumberingMng::
getCellUid(Integer level, Int64 coord_i, Int64 coord_j, Int64 coord_k)
{
  // TODO repet
  return (coord_i + coord_j * getGlobalNbCellsX(level) + coord_k * getGlobalNbCellsX(level) * getGlobalNbCellsY(level)) + getFirstCellUidLevel(level);
}

Int64 CartesianMeshNumberingMng::
getCellUid(Integer level, Int64 coord_i, Int64 coord_j)
{
  return (coord_i + coord_j * getGlobalNbCellsX(level)) + getFirstCellUidLevel(level);
}



Integer CartesianMeshNumberingMng::
getNbNode()
{
  return std::pow(m_pattern, m_mesh->dimension());
}

void CartesianMeshNumberingMng::
getNodeUids(ArrayView<Int64> uid, Integer level, Int64 coord_i, Int64 coord_j, Int64 coord_k)
{
  Int64 nb_node_x = getGlobalNbCellsX(level) + 1;
  Int64 nb_node_y = getGlobalNbCellsY(level) + 1;

  uid[0] = (coord_i + 0) + ((coord_j + 0) * nb_node_x) + ((coord_k + 0) * nb_node_x * nb_node_y);
  uid[1] = (coord_i + 1) + ((coord_j + 0) * nb_node_x) + ((coord_k + 0) * nb_node_x * nb_node_y);
  uid[2] = (coord_i + 1) + ((coord_j + 1) * nb_node_x) + ((coord_k + 0) * nb_node_x * nb_node_y);
  uid[3] = (coord_i + 0) + ((coord_j + 1) * nb_node_x) + ((coord_k + 0) * nb_node_x * nb_node_y);

  uid[4] = (coord_i + 0) + ((coord_j + 0) * nb_node_x) + ((coord_k + 1) * nb_node_x * nb_node_y);
  uid[5] = (coord_i + 1) + ((coord_j + 0) * nb_node_x) + ((coord_k + 1) * nb_node_x * nb_node_y);
  uid[6] = (coord_i + 1) + ((coord_j + 1) * nb_node_x) + ((coord_k + 1) * nb_node_x * nb_node_y);
  uid[7] = (coord_i + 0) + ((coord_j + 1) * nb_node_x) + ((coord_k + 1) * nb_node_x * nb_node_y);

  uid[0] += getFirstNodeUidLevel(level);
  uid[1] += getFirstNodeUidLevel(level);
  uid[2] += getFirstNodeUidLevel(level);
  uid[3] += getFirstNodeUidLevel(level);

  uid[4] += getFirstNodeUidLevel(level);
  uid[5] += getFirstNodeUidLevel(level);
  uid[6] += getFirstNodeUidLevel(level);
  uid[7] += getFirstNodeUidLevel(level);
}

void CartesianMeshNumberingMng::
getNodeUids(ArrayView<Int64> uid, Integer level, Int64 coord_i, Int64 coord_j)
{
  Int64 nb_node_x = getGlobalNbCellsX(level) + 1;

  uid[0] = (coord_i + 0) + ((coord_j + 0) * nb_node_x);
  uid[1] = (coord_i + 1) + ((coord_j + 0) * nb_node_x);
  uid[2] = (coord_i + 1) + ((coord_j + 1) * nb_node_x);
  uid[3] = (coord_i + 0) + ((coord_j + 1) * nb_node_x);

  uid[0] += getFirstNodeUidLevel(level);
  uid[1] += getFirstNodeUidLevel(level);
  uid[2] += getFirstNodeUidLevel(level);
  uid[3] += getFirstNodeUidLevel(level);
}

Integer CartesianMeshNumberingMng::
getNbFace()
{
  return m_pattern * m_mesh->dimension();
}

void CartesianMeshNumberingMng::
getFaceUids(ArrayView<Int64> uid, Integer level, Int64 coord_i, Int64 coord_j, Int64 coord_k)
{
  Int64 nb_cell_x = getGlobalNbCellsX(level);
  Int64 nb_cell_y = getGlobalNbCellsY(level);
  Int64 nb_cell_z = getGlobalNbCellsZ(level);

  Int64 nb_face_x = nb_cell_x + 1;
  Int64 nb_face_y = nb_cell_y + 1;
  Int64 nb_face_z = nb_cell_z + 1;

  // Numérote les faces
  // Cet algo n'est pas basé sur l'algo 2D.
  // Les UniqueIDs générés sont contigües.
  // Il est aussi possible de retrouver les UniqueIDs des faces
  // à l'aide de la position de la cellule et la taille du maillage.
  // De plus, l'ordre des UniqueIDs des faces d'une cellule est toujours le
  // même (en notation localId Arcane (cell.face(i)) : 0, 3, 1, 4, 2, 5).
  // Les UniqueIDs générés sont donc les mêmes quelque soit le découpage.
  /*
       x               z
    ┌──►          │ ┌──►
    │             │ │
   y▼12   13   14 │y▼ ┌────┬────┐
      │ 26 │ 27 │ │   │ 24 │ 25 │
      └────┴────┘ │   0    4    8
     15   16   17 │
      │ 28 │ 29 │ │   │    │    │
      └────┴────┘ │   2    5    9
   z=0            │              x=0
  - - - - - - - - - - - - - - - - - -
   z=1            │              x=1
     18   19   20 │   ┌────┬────┐
      │ 32 │ 33 │ │   │ 30 │ 31 │
      └────┴────┘ │   1    6   10
     21   22   23 │
      │ 34 │ 35 │ │   │    │    │
      └────┴────┘ │   3    7   11
                  │
  */
  // On a un cube décomposé en huit cellules (2x2x2).
  // Le schéma au-dessus représente les faces des cellules de ce cube avec
  // les uniqueIDs que l'algorithme génèrera (sans face_adder).
  // Pour cet algo, on commence par les faces "xy".
  // On énumère d'abord en x, puis en y, puis en z.
  // Une fois les faces "xy" numérotées, on fait les faces "yz".
  // Toujours le même ordre de numérotation.
  // On termine avec les faces "zx", encore dans le même ordre.
  //
  // Dans l'implémentation ci-dessous, on fait la numérotation
  // maille par maille.

  const Int64 total_face_xy = nb_face_z * nb_cell_x * nb_cell_y;
  const Int64 total_face_xy_yz = total_face_xy + nb_face_x * nb_cell_y * nb_cell_z;
  const Int64 total_face_xy_yz_zx = total_face_xy_yz + nb_face_y * nb_cell_z * nb_cell_x;

  const Int64 nb_cell_before_j = coord_j * nb_cell_x;

  uid[0] = (coord_k * nb_cell_x * nb_cell_y)
         + nb_cell_before_j
         + (coord_i);

  uid[3] = uid[0] + nb_cell_x * nb_cell_y;

  uid[1] = (coord_k * nb_face_x * nb_cell_y)
         + (coord_j * nb_face_x)
         + (coord_i) + total_face_xy;

  uid[4] = uid[1] + 1;

  uid[2] = (coord_k * nb_cell_x * nb_face_y)
         + nb_cell_before_j
         + (coord_i) + total_face_xy_yz;

  uid[5] = uid[2] + nb_cell_x;

  uid[0] += getFirstFaceUidLevel(level);
  uid[1] += getFirstFaceUidLevel(level);
  uid[2] += getFirstFaceUidLevel(level);
  uid[3] += getFirstFaceUidLevel(level);
  uid[4] += getFirstFaceUidLevel(level);
  uid[5] += getFirstFaceUidLevel(level);
}


void CartesianMeshNumberingMng::
getFaceUids(ArrayView<Int64> uid, Integer level, Int64 coord_i, Int64 coord_j)
{
  Int64 nb_cell_x = getGlobalNbCellsX(level);
  Int64 nb_face_x = nb_cell_x + 1;

  // Numérote les faces
  //  |-0--|--2-|
  // 4|   6|   8|
  //  |-5--|-7--|
  // 9|  11|  13|
  //  |-10-|-12-|
  //
  // Avec cette numérotation, HAUT < GAUCHE < BAS < DROITE
  // Mis à part les uniqueIds de la première ligne de face, tous
  // les uniqueIds sont contigües.

  // HAUT
  // - "(current_level_nb_face_x + current_level_nb_cell_x)" :
  //   le nombre de faces GAUCHE BAS DROITE au dessus.
  // - "coord_j * (current_level_nb_face_x + current_level_nb_cell_x)" :
  //   le nombre total de faces GAUCHE BAS DROITE au dessus.
  // - "coord_i * 2"
  //   on avance deux à deux sur les faces d'un même "coté".
  uid[0] = coord_i * 2 + coord_j * (nb_face_x + nb_cell_x);

  // BAS
  // Pour BAS, c'est comme HAUT mais avec un "nombre de face du dessus" en plus.
  uid[2] = uid[0] + (nb_face_x + nb_cell_x);
  // GAUCHE
  // Pour GAUCHE, c'est l'UID de BAS -1.
  uid[3] = uid[2] - 1;
  // DROITE
  // Pour DROITE, c'est l'UID de BAS +1.
  uid[1] = uid[2] + 1;

  uid[0] += getFirstFaceUidLevel(level);
  uid[1] += getFirstFaceUidLevel(level);
  uid[2] += getFirstFaceUidLevel(level);
  uid[3] += getFirstFaceUidLevel(level);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // End namespace Arcane

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
