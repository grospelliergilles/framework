//
// Created by dechaiss on 2/6/20.
//

/*-------------------------
 * Neo library
 * Polyhedral mesh test
 * sdc (C)-2020
 *
 *-------------------------
 */

#include "neo/Neo.h"
#include "gtest/gtest.h"

#ifdef  HAS_XDMF
#include "XdmfDomain.hpp"
#include "XdmfInformation.hpp"
#include "XdmfUnstructuredGrid.hpp"
#include "XdmfSystemUtils.hpp"
#include "XdmfHDF5Writer.hpp"
#include "XdmfWriter.hpp"
#include "XdmfReader.hpp"
#endif

template <typename Container>
void _printContainer(Container&& container, std::string const& name){
  std::cout << name << container.size() << std::endl;
  for (auto element : container) {
    std::cout << element << " " ;
  }
  std::cout << std::endl;
}

namespace StaticMesh {

static const std::string cell_family_name{"CellFamily"};
static const std::string face_family_name{"FaceFamily"};
static const std::string node_family_name{"NodeFamily"};

void addItems(Neo::Mesh& mesh, Neo::Family& family, std::vector<Neo::utils::Int64> const& uids, Neo::AddedItemRange& added_item_range)
{
  auto& added_items = added_item_range.new_items;
  // Add items
  mesh.addAlgorithm(  Neo::OutProperty{family,family.lidPropName()},
                      [&family,&uids,&added_items](Neo::ItemLidsProperty & lids_property){
                        std::cout << "Algorithm: create items in family " << family.m_name << std::endl;
                        added_items = lids_property.append(uids);
                        lids_property.debugPrint();
                        std::cout << "Inserted item range : " << added_items;
                      });
  // register their uids
  auto uid_property_name = family.m_name+"_uids";
  mesh.addAlgorithm(
      Neo::InProperty{family,family.lidPropName()},
      Neo::OutProperty{family, uid_property_name},
      [&family,&uids,&added_items](Neo::ItemLidsProperty const& item_lids_property,
                                   Neo::PropertyT<Neo::utils::Int64>& item_uids_property){
        std::cout << "Algorithm: register item uids for family " << family.m_name << std::endl;
        if (item_uids_property.isInitializableFrom(added_items))
          item_uids_property.init(added_items,std::move(uids)); // init can steal the input values
        else
          item_uids_property.append(added_items, uids);
        item_uids_property.debugPrint();
      });// need to add a property check for existing uid
  }

  // todo same interface with nb_connected_item_per_items as an array
  void addConnectivity(Neo::Mesh &mesh, Neo::Family &source_family,
                       Neo::ItemRange &source_items,
                       Neo::Family& target_family,
                       std::vector<size_t>&& nb_connected_item_per_items,
                       std::vector<Neo::utils::Int64> const& connected_item_uids) {
    // add connectivity property if doesn't exist
    std::string connectivity_name =
        source_family.m_name + "to" + target_family.m_name + "_connectivity";
    source_family.addArrayProperty<Neo::utils::Int32>(connectivity_name);
    mesh.addAlgorithm(
        Neo::InProperty{source_family, source_family.lidPropName()},
        Neo::InProperty{target_family, target_family.lidPropName()},
        Neo::OutProperty{source_family, source_family.m_name + "to" +
                                            target_family.m_name +
                                            "_connectivity"},
        [&connected_item_uids, nb_connected_item_per_items, &source_items,
         &source_family, &target_family](
            Neo::ItemLidsProperty const &source_family_lids_property,
            Neo::ItemLidsProperty const &target_family_lids_property,
            Neo::ArrayProperty<Neo::utils::Int32> &source2target) {
          std::cout << "Algorithm: register connectivity between "
                    << source_family.m_name << "  and  " << target_family.m_name
                    << std::endl;
          auto connected_item_lids =
              target_family_lids_property[connected_item_uids];
          if (source2target.isInitializableFrom(source_items)) {
            source2target.resize(std::move(nb_connected_item_per_items));
            source2target.init(source_items, std::move(connected_item_lids));
          } else {
            source2target.append(source_items, connected_item_lids,
                                 nb_connected_item_per_items);
          }
          source2target.debugPrint();
        });
  }

  Neo::ArrayProperty<Neo::utils::Int32> const &
  getConnectivity(Neo::Mesh const &mesh, Neo::Family const& source_family,
                  Neo::Family const &target_family)
  {
    return source_family.getConcreteProperty <
           Neo::ArrayProperty<Neo::utils::Int32>>(source_family.m_name + "to" +
                                                 target_family.m_name+"_connectivity");
  }

  // todo : define 2 signatures to indicate eventual memory stealing...?
  void setNodeCoords(Neo::Mesh& mesh, Neo::Family& node_family, Neo::AddedItemRange& added_node_range, std::vector<Neo::utils::Real3>& node_coords){
    node_family.addProperty<Neo::utils::Real3>(std::string("node_coords"));
    auto& added_nodes = added_node_range.new_items;
    mesh.addAlgorithm(
        Neo::InProperty{node_family,node_family.lidPropName()},
        Neo::OutProperty{node_family,"node_coords"},
        [&node_coords,&added_nodes](Neo::ItemLidsProperty const& node_lids_property,
                                    Neo::PropertyT<Neo::utils::Real3> & node_coords_property){
          std::cout << "Algorithm: register node coords" << std::endl;
          if (node_coords_property.isInitializableFrom(added_nodes))  node_coords_property.init(added_nodes,std::move(node_coords)); // init can steal the input values
          else node_coords_property.append(added_nodes, node_coords);
          node_coords_property.debugPrint();
        });
  }

  // todo define also a const method
  Neo::utils::ArrayView<Neo::utils::Real3> getNodeCoords(Neo::Mesh const& mesh, Neo::Family& node_family)
  {
    auto& node_coords = node_family.getConcreteProperty<Neo::PropertyT<Neo::utils::Real3>>("node_coords");
    return node_coords.values(); //
  }

  void addConnectivity(Neo::Mesh &mesh, Neo::Family &source_family,
                       Neo::AddedItemRange &source_items,
                       Neo::Family& target_family,
                       std::vector<size_t>&& nb_connected_item_per_items,
                       std::vector<Neo::utils::Int64> const& connected_item_uids)
  {
    addConnectivity(mesh, source_family, source_items.new_items, target_family,
                    std::move(nb_connected_item_per_items), connected_item_uids);
  }

  Neo::ArrayProperty<Neo::utils::Int32> const& faces(Neo::Mesh const& mesh, Neo::Family const& source_family)
  {
    auto &face_family = mesh.getFamily(Neo::ItemKind::IK_Face, face_family_name);
    return getConnectivity(mesh, source_family, face_family);
  }

  Neo::ArrayProperty<Neo::utils::Int32> const& nodes(Neo::Mesh const& mesh, Neo::Family const& source_family)
  {
    auto& node_family = mesh.getFamily(Neo::ItemKind::IK_Node,node_family_name);
    return getConnectivity(mesh, source_family, node_family);
  }
}

namespace PolyhedralMeshTest {

  auto &addCellFamily(Neo::Mesh &mesh, std::string family_name) {
    auto &cell_family =
        mesh.addFamily(Neo::ItemKind::IK_Cell, std::move(family_name));
    cell_family.addProperty<Neo::utils::Int64>(family_name + "_uids");
    return cell_family;
  }

  auto &addNodeFamily(Neo::Mesh &mesh, std::string family_name) {
    auto &node_family =
        mesh.addFamily(Neo::ItemKind::IK_Node, std::move(family_name));
    node_family.addProperty<Neo::utils::Int64>(family_name + "_uids");
    return node_family;
  }

  auto &addFaceFamily(Neo::Mesh &mesh, std::string family_name) {
    auto &face_family =
        mesh.addFamily(Neo::ItemKind::IK_Face, std::move(family_name));
    face_family.addProperty<Neo::utils::Int64>(family_name + "_uids");
    return face_family;
  }

  void _createMesh(Neo::Mesh &mesh,
                   std::vector<Neo::utils::Int64> const &node_uids,
                   std::vector<Neo::utils::Int64> const &cell_uids,
                   std::vector<Neo::utils::Int64> const &face_uids,
                   std::vector<Neo::utils::Real3>& node_coords, // not const since they can be moved
                   std::vector<Neo::utils::Int64>& cell_nodes,
                   std::vector<Neo::utils::Int64>& cell_faces,
                   std::vector<Neo::utils::Int64>& face_nodes,
                   std::vector<size_t>&& nb_node_per_cells,
                   std::vector<size_t>&& nb_node_per_faces,
                   std::vector<size_t>&& nb_face_per_cells) {
    auto &cell_family = addCellFamily(mesh, StaticMesh::cell_family_name);
    auto &node_family = addNodeFamily(mesh, StaticMesh::node_family_name);
    auto &face_family = addFaceFamily(mesh, StaticMesh::face_family_name);
    mesh.beginUpdate();
    auto added_cells = Neo::AddedItemRange{};
    auto added_nodes = Neo::AddedItemRange{};
    auto added_faces = Neo::AddedItemRange{};
    StaticMesh::addItems(mesh, cell_family, cell_uids, added_cells);
    StaticMesh::addItems(mesh, node_family, node_uids, added_nodes);
    StaticMesh::addItems(mesh, face_family, face_uids, added_faces);
    StaticMesh::setNodeCoords(mesh, node_family, added_nodes, node_coords);
    StaticMesh::addConnectivity(mesh, cell_family, added_cells, node_family,
                                std::move(nb_node_per_cells), cell_nodes);
    StaticMesh::addConnectivity(mesh, face_family, added_faces, node_family,
                                std::move(nb_node_per_faces), face_nodes);
    StaticMesh::addConnectivity(mesh, cell_family, added_cells, face_family,
                                std::move(nb_face_per_cells), cell_faces);
    auto valid_mesh_state = mesh.endUpdate();
    auto &new_cells = added_cells.get(valid_mesh_state);
    auto &new_nodes = added_nodes.get(valid_mesh_state);
    auto &new_faces = added_faces.get(valid_mesh_state);
    std::cout << "Added cells range after endUpdate: " << new_cells;
    std::cout << "Added nodes range after endUpdate: " << new_nodes;
    std::cout << "Added faces range after endUpdate: " << new_faces;
  }

  void addCells(Neo::Mesh &mesh) {
    std::vector<Neo::utils::Int64> node_uids{0, 1, 2, 3, 4, 5};
    std::vector<Neo::utils::Int64> cell_uids{0};
    std::vector<Neo::utils::Int64> face_uids{0, 1, 2, 3, 4, 5, 6, 7};

    std::vector<Neo::utils::Real3> node_coords{
        {-1, -1, 0}, {-1, 1, 0}, {1, 1, 0}, {1, -1, 0}, {0, 0, 1}, {0, 0, -1}};

    std::vector<Neo::utils::Int64> cell_nodes{0, 1, 2, 3, 4, 5};
    std::vector<Neo::utils::Int64> cell_faces{0, 1, 2, 3, 4, 5, 6, 7};

    std::vector<Neo::utils::Int64> face_nodes{0, 1, 4, 0, 1, 5, 1, 2, 4, 1, 2, 5,
                                              2, 3, 4, 2, 3, 5, 3, 0, 4, 3, 0, 5};

    auto nb_node_per_cell = 6;
    auto nb_node_per_face = 3;
    auto nb_face_per_cell = 8;

    _createMesh(mesh, node_uids, cell_uids, face_uids, node_coords, cell_nodes,
                cell_faces, face_nodes,
                std::vector<size_t>(cell_uids.size(),nb_node_per_cell),
                std::vector<size_t>(face_uids.size(),nb_node_per_face),
                std::vector<size_t>(cell_uids.size(),nb_face_per_cell));
  }
}

namespace XdmfTest {
  void exportMesh(Neo::Mesh const& mesh, std::string const& file_name)
  {
    auto domain = XdmfDomain::New();
    auto domain_info = XdmfInformation::New("Domain", " For polyhedral data from Neo");
    domain->insert(domain_info);
    // Needed ?
    auto xdmf_grid = XdmfUnstructuredGrid::New();
    auto xdmf_geom = XdmfGeometry::New();
    xdmf_geom->setType(XdmfGeometryType::XYZ());
    auto node_coords = StaticMesh::getNodeCoords(mesh,mesh.getFamily(Neo::ItemKind::IK_Node,StaticMesh::node_family_name));
    xdmf_geom->insert(0,(double*)node_coords.begin(),node_coords.size()*3,1,1);
    xdmf_grid->setGeometry(xdmf_geom);
    auto xdmf_topo = XdmfTopology::New();
    xdmf_topo->setType(XdmfTopologyType::Polyhedron());
    std::vector<int> cell_data;
    auto& cell_family = mesh.getFamily(Neo::ItemKind::IK_Cell,StaticMesh::cell_family_name);
    auto& face_family = mesh.getFamily(Neo::ItemKind::IK_Face,StaticMesh::face_family_name);
    cell_data.reserve(cell_family.nbElements()*4); // 4 faces by cell approx
    auto& cell_to_faces = StaticMesh::faces(mesh, cell_family);
    auto& face_to_nodes = StaticMesh::nodes(mesh, face_family);
    for (auto cell : cell_family.all()) {
      auto cell_faces = cell_to_faces[cell];
      cell_data.push_back(cell_faces.size());
      for (auto face : cell_faces) {
        auto face_nodes = face_to_nodes[face];
        cell_data.push_back(face_nodes.size());
        cell_data.insert(cell_data.end(),face_nodes.begin(),face_nodes.end());
      }
    }
    xdmf_topo->insert(0, cell_data.data(),cell_data.size(), 1, 1);
    xdmf_grid->setTopology(xdmf_topo);
    domain->insert(xdmf_grid);
//  auto heavydata_writer = XdmfHDF5Writer::New(file_name);
//  auto writer = XdmfWriter::New(file_name, heavydata_writer);
    auto writer = XdmfWriter::New(file_name);
//  domain->accept(heavydata_writer);
    writer->setLightDataLimit(1000);
    domain->accept(writer);
  }
}

TEST(PolyhedralTest,CreateMesh1)
{
  auto mesh = Neo::Mesh{"PolyhedralMesh"};
  PolyhedralMeshTest::addCells(mesh);
}


#ifdef HAS_XDMF
TEST(PolyhedralTest,CreateXdmfMesh)
{
  auto mesh = Neo::Mesh{"PolyhedralMesh"};
  PolyhedralMeshTest::addCells(mesh);
  auto exported_mesh{"test_output.xmf"};
  XdmfTest::exportMesh(mesh,exported_mesh);
  // todo reimport to check
  auto reader = XdmfReader::New();
  auto exported_primaryDomain = shared_dynamic_cast<XdmfDomain>(reader->read(exported_mesh));
  auto ref_primaryDomain = shared_dynamic_cast<XdmfDomain>(reader->read("../test/meshes/example_cell.xmf"));
  auto exported_topology_str = exported_primaryDomain->getUnstructuredGrid("Grid")->getTopology()->getValuesString();
  auto ref_topology_str = ref_primaryDomain->getUnstructuredGrid("Octahedron")->getTopology()->getValuesString();
  auto exported_geometry_str = exported_primaryDomain->getUnstructuredGrid("Grid")->getGeometry()->getValuesString();
  auto ref_geometry_str = ref_primaryDomain->getUnstructuredGrid("Octahedron")->getGeometry()->getValuesString();
  std::cout << "original topology " << ref_topology_str << std::endl;
  std::cout << "exported topology " << exported_topology_str << std::endl;
  std::cout << "original geometry " << ref_geometry_str<< std::endl;
  std::cout << "exported geometry " << exported_geometry_str << std::endl;
  EXPECT_EQ(std::string{ref_topology_str},std::string{exported_topology_str}); // comparer avec std::equal
  EXPECT_EQ(std::string{ref_geometry_str},std::string{exported_geometry_str}); // comparer avec std::equal

}

TEST(PolyhedralTest,ImportXdmfMesh)
{
  auto reader = XdmfReader::New();
  auto primaryDomain = shared_dynamic_cast<XdmfDomain>(reader->read("../test/meshes/example_mesh.xmf"));
  auto grid = primaryDomain->getUnstructuredGrid("Polyhedra");
  auto geometry = grid->getGeometry();
  geometry->read();
  EXPECT_EQ(geometry->getType()->getName(), XdmfGeometryType::XYZ()->getName());
  std::vector<Neo::utils::Real3> node_coords(geometry->getNumberPoints(),{-1e6,-1e6,-1e6});
  geometry->getValues(0,(double*)node_coords.data(),geometry->getNumberPoints()*3,1,1);
  auto topology = grid->getTopology();
  topology->read();
  // Read only polyhedrons
  EXPECT_EQ(XdmfTopologyType::Polyhedron()->getName(),topology->getType()->getName());
  std::vector<Neo::utils::Int32> cell_data(topology->getSize(),-1);
  topology->getValues(0,cell_data.data(),topology->getSize());
  std::vector<Neo::utils::Int64> cell_uids;
  std::vector<Neo::utils::Int64> face_uids;
  std::vector<Neo::utils::Int64> face_nodes;
  std::set<Neo::utils::Int64> node_uids_set;
  std::vector<Neo::utils::Int64> node_uids;
  std::set<Neo::utils::Int32> current_cell_nodes;
  std::vector<Neo::utils::Int64> cell_nodes;
  std::vector<Neo::utils::Int64> cell_faces;
  std::vector<size_t> nb_node_per_cells;
  std::vector<size_t> nb_face_per_cells;
  using FaceNodes = std::set<int>;
  using FaceUid = Neo::utils::Int64;
  using FaceInfo = std::pair<FaceNodes,FaceUid>;
  auto face_info_comp = [](FaceInfo const& face_info1, FaceInfo const& face_info2){
    return face_info1.first < face_info2.first;
  };
  std::set<FaceInfo, decltype(face_info_comp)> face_nodes_set(face_info_comp);
  face_nodes.reserve(topology->getSize());
  std::vector<size_t> nb_node_per_faces;
  auto cell_index = 0;
  auto face_uid = 0;
  for (auto cell_data_index = 0; cell_data_index< cell_data.size();){
    cell_uids.push_back(cell_index++);
    auto cell_nb_face  = cell_data[cell_data_index++];
    nb_face_per_cells.push_back(cell_nb_face);
    for (auto face_index = 0; face_index< cell_nb_face;++face_index){
      std::size_t face_nb_node = cell_data[cell_data_index++];
      auto current_face_nodes = Neo::utils::ArrayView<Neo::utils::Int32>{face_nb_node,&cell_data[cell_data_index]};
      auto [face_info, is_new_face] = face_nodes_set.emplace(FaceNodes{current_face_nodes.begin(),
                             current_face_nodes.end()},face_uid);
      if (!is_new_face) std::cout << "Face not inserted " << face_uid << std::endl;
      if (is_new_face) {
        face_nodes.insert(face_nodes.end(),current_face_nodes.begin(), current_face_nodes.end());
        nb_node_per_faces.push_back(face_nb_node);
      }
      cell_faces.push_back(face_info->second);
      if (is_new_face) face_uids.push_back(face_uid++);
      current_cell_nodes.insert(current_face_nodes.begin(), current_face_nodes.end());
      cell_data_index += face_nb_node;
    }
    cell_nodes.insert(cell_nodes.end(), current_cell_nodes.begin(),current_cell_nodes.end());
    nb_node_per_cells.push_back(current_cell_nodes.size());
    node_uids_set.insert(current_cell_nodes.begin(),current_cell_nodes.end());
    current_cell_nodes.clear();
  }
  node_uids.insert(node_uids.end(),std::begin(node_uids_set), std::end(node_uids_set));
  _printContainer(face_nodes, "face nodes ");
  _printContainer(face_uids, "face uids ");
  _printContainer(nb_node_per_faces, "nb node per face ");
  _printContainer(node_uids_set, "node uids ");
  _printContainer(cell_nodes, "cell nodes ");
  _printContainer(nb_node_per_cells, "nb node per cell ");
  _printContainer(cell_faces, "cell faces ");
  _printContainer(nb_face_per_cells, "nb face per cell ");
  // local checks
  std::vector<Neo::utils::Int64> cell_uids_ref = {0, 1, 2};
  EXPECT_TRUE(std::equal(cell_uids.begin(),cell_uids.end(),cell_uids_ref.begin(),cell_uids_ref.end()));
  EXPECT_EQ(27,face_uids.size());
  EXPECT_EQ(geometry->getNumberPoints(), node_uids_set.size());
  // import mesh in Neo data structure
  auto mesh = Neo::Mesh{"'ImportedMesh"};
  PolyhedralMeshTest::_createMesh(mesh, node_uids,cell_uids,face_uids,node_coords,cell_nodes,cell_faces,face_nodes,
      std::move(nb_node_per_cells),std::move(nb_node_per_faces),std::move(nb_face_per_cells));
  std::string imported_mesh{"imported_mesh.xmf"};
  XdmfTest::exportMesh(mesh,imported_mesh);
  // Compare with original mesh
  auto created_primaryDomain = shared_dynamic_cast<XdmfDomain>(reader->read(imported_mesh));
  std::cout << "original topology " << topology->getValuesString().c_str()<< std::endl;
  std::cout << "created topology " << created_primaryDomain->getUnstructuredGrid("Grid")->getTopology()->getValuesString().c_str()<< std::endl;
  std::cout << "original geometry " << geometry->getValuesString().c_str()<< std::endl;
  std::cout << "created geometry " << created_primaryDomain->getUnstructuredGrid("Grid")->getGeometry()->getValuesString().c_str()<< std::endl;
  EXPECT_EQ(std::string{geometry->getValuesString().c_str()},std::string{created_primaryDomain->getUnstructuredGrid("Grid")->getGeometry()->getValuesString().c_str()}); // comparer avec std::equal
  EXPECT_EQ(std::string{topology->getValuesString().c_str()},std::string{created_primaryDomain->getUnstructuredGrid("Grid")->getTopology()->getValuesString().c_str()}); // comparer avec std::equal
}

#endif