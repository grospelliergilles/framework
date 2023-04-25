﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* CustomMeshTestModule.cc                          C) 2000-2021             */
/*                                                                           */
/* Test Module for custom mesh                                               */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <numeric>

#include "arcane/ITimeLoopMng.h"
#include "arcane/IMesh.h"
#include "arcane/MeshHandle.h"
#include "arcane/IMeshMng.h"
#include "arcane/mesh/PolyhedralMesh.h"
#include "arcane/utils/ValueChecker.h"

#include "arcane/ItemGroup.h"
#include "CustomMeshTest_axl.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace ArcaneTest::CustomMesh
{
using namespace Arcane;
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

class CustomMeshTestModule : public ArcaneCustomMeshTestObject
{
 public:

  explicit CustomMeshTestModule(const ModuleBuildInfo& sbi)
  : ArcaneCustomMeshTestObject(sbi)
  {}

 public:

  void init()
  {
    info() << "-- INIT CUSTOM MESH MODULE";
    auto mesh_handle = subDomain()->defaultMeshHandle();
    if (mesh_handle.hasMesh()) {
      info() << "-- MESH NAME: " << mesh()->name();
      _testDimensions(mesh());
      _testCoordinates(mesh());
      _testEnumerationAndConnectivities(mesh());
      _testVariables(mesh());
      _testGroups(mesh());
    }
    else
      info() << "No Mesh";

    subDomain()->timeLoopMng()->stopComputeLoop(true);
  }

 private:

  void _testEnumerationAndConnectivities(IMesh* mesh);
  void _testVariables(IMesh* mesh);
  void _testGroups(IMesh* mesh);
  void _testDimensions(IMesh* mesh);
  void _testCoordinates(IMesh* mesh);
  void _buildGroup(IItemFamily* family, String const& group_name);
  template <typename VariableRefType>
  void _checkVariable(VariableRefType variable, ItemGroup item_group);
  template <typename VariableArrayRefType>
  void _checkArrayVariable(VariableArrayRefType variable, ItemGroup item_group);
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CustomMeshTestModule::
_testEnumerationAndConnectivities(IMesh* mesh)
{
  info() << "- Polyhedral mesh test -";
  info() << "- Mesh dimension " << mesh->dimension();
  info() << "- Mesh nb cells  " << mesh->nbItem(IK_Cell) << " or " << mesh->nbCell();
  info() << "- Mesh nb faces  " << mesh->nbItem(IK_Face) << " or " << mesh->nbFace();
  info() << "- Mesh nb edges  " << mesh->nbItem(IK_Edge) << " or " << mesh->nbEdge();
  info() << "- Mesh nb nodes  " << mesh->nbItem(IK_Node) << " or " << mesh->nbNode();
  info() << "Cell family " << mesh->cellFamily()->name();
  info() << "Node family " << mesh->nodeFamily()->name();
  auto all_cells = mesh->allCells();
  // ALL CELLS
  ENUMERATE_CELL (icell, all_cells) {
    info() << "cell with index " << icell.index();
    info() << "cell with lid " << icell.localId();
    info() << "cell with uid " << icell->uniqueId().asInt64();
    info() << "cell number of nodes " << icell->nodes().size();
    info() << "cell number of faces " << icell->faces().size();
    info() << "cell number of edges " << icell->edges().size();
    ENUMERATE_NODE (inode, icell->nodes()) {
      info() << "cell node " << inode.index() << " lid " << inode.localId() << " uid " << inode->uniqueId().asInt64();
    }
    ENUMERATE_FACE (iface, icell->faces()) {
      info() << "cell face " << iface.index() << " lid " << iface.localId() << " uid " << iface->uniqueId().asInt64();
    }
    ENUMERATE_EDGE (iedge, icell->edges()) {
      info() << "cell edge " << iedge.index() << " lid " << iedge.localId() << " uid " << iedge->uniqueId().asInt64();
    }
  }
  // ALL FACES
  ENUMERATE_FACE (iface, mesh->allFaces()) {
    info() << "face with index " << iface.index();
    info() << "face with lid " << iface.localId();
    info() << "face with uid " << iface->uniqueId().asInt64();
    info() << "face number of nodes " << iface->nodes().size();
    info() << "face number of cells " << iface->cells().size();
    info() << "face number of edges " << iface->edges().size();
    ENUMERATE_NODE (inode, iface->nodes()) {
      info() << "face node " << inode.index() << " lid " << inode.localId() << " uid " << inode->uniqueId().asInt64();
    }
    ENUMERATE_CELL (icell, iface->cells()) {
      info() << "face cell " << icell.index() << " lid " << icell.localId() << " uid " << icell->uniqueId().asInt64();
    }
    ENUMERATE_EDGE (iedge, iface->edges()) {
      info() << "face edge " << iedge.index() << " lid " << iedge.localId() << " uid " << iedge->uniqueId().asInt64();
    }
  }
  // ALL NODES
  ENUMERATE_NODE (inode, mesh->allNodes()) {
    info() << "node with index " << inode.index();
    info() << "node with lid " << inode.localId();
    info() << "node with uid " << inode->uniqueId().asInt64();
    info() << "node number of faces " << inode->faces().size();
    info() << "node number of cells " << inode->cells().size();
    info() << "node number of edges " << inode->edges().size();
    ENUMERATE_FACE (iface, inode->faces()) {
      info() << "node face " << iface.index() << " lid " << iface.localId() << " uid " << iface->uniqueId().asInt64();
    }
    ENUMERATE_CELL (icell, inode->cells()) {
      info() << "node cell " << icell.index() << " lid " << icell.localId() << " uid " << icell->uniqueId().asInt64();
    }
    ENUMERATE_EDGE (iedge, inode->edges()) {
      info() << "face edge " << iedge.index() << " lid " << iedge.localId() << " uid " << iedge->uniqueId().asInt64();
    }
  }
  // ALL EDGES
  ENUMERATE_EDGE (iedge, mesh->allEdges()) {
    info() << "edge with index " << iedge.index();
    info() << "edge with lid " << iedge.localId();
    info() << "edge with uid " << iedge->uniqueId().asInt64();
    info() << "edge number of faces " << iedge->faces().size();
    info() << "edge number of cells " << iedge->cells().size();
    info() << "edge number of nodes " << iedge->nodes().size();
    ENUMERATE_FACE (iface, iedge->faces()) {
      info() << "edge face " << iface.index() << " lid " << iface.localId() << " uid " << iface->uniqueId().asInt64();
    }
    ENUMERATE_CELL (icell, iedge->cells()) {
      info() << "edge cell " << icell.index() << " lid " << icell.localId() << " uid " << icell->uniqueId().asInt64();
    }
    ENUMERATE_NODE (inode, iedge->nodes()) {
      info() << "edge node " << inode.index() << " lid " << inode.localId() << " uid " << inode->uniqueId().asInt64();
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CustomMeshTestModule::_testVariables(IMesh* mesh)
{
  // test variables
  info() << " -- test variables -- ";
  // cell variable
  m_cell_variable.fill(1);
  _checkVariable(m_cell_variable, mesh->allCells());
  // node variable
  m_node_variable.fill(1);
  _checkVariable(m_node_variable, mesh->allNodes());
  // face variable
  m_face_variable.fill(1);
  _checkVariable(m_face_variable, mesh->allFaces());
  // edge variable
  m_edge_variable.fill(1);
  _checkVariable(m_edge_variable, mesh->allEdges());
  // Check variables defined in mesh file
  // Cell variables
  for (const auto& variable_name : options()->getCheckCellVariableReal()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableCellReal var{ VariableBuildInfo(mesh, variable_name) };
    _checkVariable(var, mesh->allCells());
  }
  for (const auto& variable_name : options()->getCheckCellVariableInteger()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableCellInteger var{ VariableBuildInfo(mesh, variable_name) };
    _checkVariable(var, mesh->allCells());
  }
  for (const auto& variable_name : options()->getCheckCellVariableArrayInteger()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableCellArrayInteger var{ VariableBuildInfo(mesh, variable_name) };
    _checkArrayVariable(var, mesh->allCells());
  }
  for (const auto& variable_name : options()->getCheckCellVariableArrayReal()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableCellArrayReal var{ VariableBuildInfo(mesh, variable_name) };
    _checkArrayVariable(var, mesh->allCells());
  }
  // Node variables
  for (const auto& variable_name : options()->getCheckNodeVariableReal()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableNodeReal var{ VariableBuildInfo(mesh, variable_name) };
    _checkVariable(var, mesh->allNodes());
  }
  for (const auto& variable_name : options()->getCheckNodeVariableInteger()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableNodeInteger var{ VariableBuildInfo(mesh, variable_name) };
    _checkVariable(var, mesh->allNodes());
  }
  for (const auto& variable_name : options()->getCheckNodeVariableArrayInteger()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableNodeArrayInteger var{ VariableBuildInfo(mesh, variable_name) };
    _checkArrayVariable(var, mesh->allNodes());
  }
  for (const auto& variable_name : options()->getCheckNodeVariableArrayReal()) {
    if (!Arcane::AbstractModule::subDomain()->variableMng()->findMeshVariable(mesh, variable_name))
      ARCANE_FATAL("Cannot find mesh array variable {0}", variable_name);
    VariableNodeArrayReal var{ VariableBuildInfo(mesh, variable_name) };
    _checkArrayVariable(var, mesh->allNodes());
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CustomMeshTestModule::
_testGroups(IMesh* mesh)
{
  // AllItems groups
  ARCANE_ASSERT((!mesh->findGroup("AllCells").null()), ("Group AllCells has not been created"));
  ARCANE_ASSERT((!mesh->findGroup("AllNodes").null()), ("Group AllNodes has not been created"));
  ARCANE_ASSERT((!mesh->findGroup("AllFaces").null()), ("Group AllFaces has not been created"));
  ARCANE_ASSERT((!mesh->findGroup("AllEdges").null()), ("Group AllEdges has not been created"));
  // Cell group
  String group_name = "my_cell_group";
  _buildGroup(mesh->cellFamily(), group_name);
  ARCANE_ASSERT((!mesh->findGroup(group_name).null()), ("Group my_cell_group has not been created"));
  PartialVariableCellInt32 partial_cell_var({ mesh, "partial_cell_variable", mesh->cellFamily()->name(), group_name });
  partial_cell_var.fill(1);
  _checkVariable(partial_cell_var, partial_cell_var.itemGroup());
  // Node group
  group_name = "my_node_group";
  _buildGroup(mesh->nodeFamily(), group_name);
  ARCANE_ASSERT((!mesh->findGroup(group_name).null()), ("Group my_node_group has not been created"));
  PartialVariableNodeInt32 partial_node_var({ mesh, "partial_node_variable", mesh->nodeFamily()->name(), group_name });
  partial_node_var.fill(1);
  _checkVariable(partial_node_var, partial_node_var.itemGroup());
  // Face group
  group_name = "my_face_group";
  _buildGroup(mesh->faceFamily(), group_name);
  ARCANE_ASSERT((!mesh->findGroup(group_name).null()), ("Group my_face_group has not been created"));
  PartialVariableFaceInt32 partial_face_var({ mesh, "partial_face_variable", mesh->faceFamily()->name(), group_name });
  partial_face_var.fill(1);
  _checkVariable(partial_face_var, partial_face_var.itemGroup());
  // Edge group
  group_name = "my_edge_group";
  _buildGroup(mesh->edgeFamily(), group_name);
  ARCANE_ASSERT((!mesh->findGroup(group_name).null()), ("Group my_edge_group has not been created"));
  PartialVariableEdgeInt32 partial_edge_var({ mesh, "partial_edge_variable", mesh->edgeFamily()->name(), group_name });
  partial_edge_var.fill(1);
  _checkVariable(partial_edge_var, partial_edge_var.itemGroup());

  for (const auto& group_infos : options()->checkGroup()) {
    auto group = mesh->findGroup(group_infos->getName());
    if (group.null())
      ARCANE_FATAL("Could not find group {0}", group_infos->getName());
    ValueChecker vc{ A_FUNCINFO };
    vc.areEqual(group.size(), group_infos->getSize(), "check group size");
  }
  ValueChecker vc{ A_FUNCINFO };
  auto nb_group = 8 + options()->checkGroup().size();
  vc.areEqual(nb_group, mesh->groups().count(), "check number of groups in the mesh");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CustomMeshTestModule::
_testDimensions(IMesh* mesh)
{
  auto mesh_size = options()->meshSize();
  if (mesh_size.empty())
    return;
  ValueChecker vc(A_FUNCINFO);
  vc.areEqual(mesh->nbCell(), mesh_size[0]->getNbCells(), "check number of cells");
  vc.areEqual(mesh->nbFace(), mesh_size[0]->getNbFaces(), "check number of faces");
  vc.areEqual(mesh->nbEdge(), mesh_size[0]->getNbEdges(), "check number of edges");
  vc.areEqual(mesh->nbNode(), mesh_size[0]->getNbNodes(), "check number of nodes");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CustomMeshTestModule::
_testCoordinates(Arcane::IMesh* mesh)
{
  if (options()->meshCoordinates.size() == 1) {
    if (options()->meshCoordinates[0].doCheck()) {
      auto node_coords = mesh->toPrimaryMesh()->nodesCoordinates();
      auto node_coords_ref = options()->meshCoordinates[0].coords();
      ValueChecker vc{ A_FUNCINFO };
      ENUMERATE_NODE (inode, allNodes()) {
        vc.areEqual(node_coords[inode], node_coords_ref[0]->value[inode.index()], "check coords values");
        info() << " node coords  " << node_coords[inode];
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void CustomMeshTestModule::
_buildGroup(IItemFamily* family, String const& group_name)
{
  auto group = family->findGroup(group_name, true);
  Int32UniqueArray item_lids;
  item_lids.reserve(family->nbItem());
  ENUMERATE_ITEM (iitem, family->allItems()) {
    if (iitem.localId() % 2 == 0)
      item_lids.add(iitem.localId());
  }
  group.addItems(item_lids);
  info() << itemKindName(family->itemKind()) << " group size " << group.size();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <typename VariableRefType>
void CustomMeshTestModule::
_checkVariable(VariableRefType variable_ref, ItemGroup item_group)
{
  auto variable_sum = 0.;
  using ItemType = typename VariableRefType::ItemType;
  ENUMERATE_ (ItemType, iitem, item_group) {
    info() << variable_ref.name() << " at item " << iitem.localId() << " " << variable_ref[iitem];
    variable_sum += variable_ref[iitem];
  }
  if (variable_sum != item_group.size()) {
    fatal() << "Error on variable " << variable_ref.name();
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template <typename VariableArrayRefType>
void CustomMeshTestModule::
_checkArrayVariable(VariableArrayRefType variable_ref, ItemGroup item_group)
{
  auto variable_sum = 0.;
  using ItemType = typename VariableArrayRefType::ItemType;
  auto array_size = variable_ref.arraySize();
  ENUMERATE_ (ItemType, iitem, item_group) {
    for (auto value : variable_ref[iitem]) {
      variable_sum += value;
    }
    info() << variable_ref.name() << " at item " << iitem.localId() << variable_ref[iitem];
  }
  ValueChecker vc{ A_FUNCINFO };
  if (array_size == 0)
    ARCANE_FATAL("Array variable {0} array size is zero");
  std::vector<int> ref_sum(array_size);
  std::iota(ref_sum.begin(), ref_sum.end(), 1.);
  vc.areEqual(variable_sum, item_group.size() * std::accumulate(ref_sum.begin(), ref_sum.end(), 0.), "check array variable values");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_REGISTER_MODULE_CUSTOMMESHTEST(CustomMeshTestModule);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // End namespace ArcaneTest::CustomMesh

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/