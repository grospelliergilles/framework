// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArcanePybind11ArrayData.cpp                             (C) 2000-2026 */
/*                                                                           */
/* Implementation of pybind11 bindings for ArrayDataT.                       */
/*---------------------------------------------------------------------------*/

#include "arcane/python/ArrayDataPybind11.h"

#include <sstream>

namespace py = pybind11;

namespace Arcane::Pybind11
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Implementation of ArrayDataT wrapper
/*---------------------------------------------------------------------------*/

template<typename DataType>
void register_ArrayDataT_impl(py::module_& m, const std::string& type_name)
{
  using ThisClass = ArrayDataT<DataType>;
  using DataInterfaceType = IArrayDataT<DataType>;

  std::string class_name = type_name + "ArrayData";

  py::class_<ThisClass, DataInterfaceType, IData, std::unique_ptr<ThisClass, py::deleter_base>>(m, class_name.c_str())
    .def(py::init<ITraceMng*>(), py::arg("trace"), "Constructor with trace")
    .def(py::init<const DataStorageBuildInfo&>(), py::arg("dsbi"), "Constructor with build info")
    .def(py::init<const ThisClass&>(), py::arg("rhs"), "Copy constructor")
    .def("__copy__", [](const ThisClass& self) { return new ThisClass(self); })
    .def("__deepcopy__", [](const ThisClass& self, py::object) { return new ThisClass(self); })
    .def("dimension", &ThisClass::dimension, "Return dimension (always 1)")
    .def("multiTag", &ThisClass::multiTag, "Return multi-tag (always 0)")
    .def("dataType", &ThisClass::dataType, "Return data type")
    .def("view", py::overload_cast<>(&ThisClass::view), "Get mutable view on data")
    .def("view", py::overload_cast<const ThisClass*>(&ThisClass::view, py::const), "Get const view on data")
    .def("value", py::overload_cast<>(&ThisClass::value), "Get mutable reference to data array")
    .def("value", py::overload_cast<const ThisClass*>(&ThisClass::value, py::const), "Get const reference to data array")
    .def("resize", &ThisClass::resize, py::arg("new_size"), "Resize the data array")
    .def("cloneRef", py::overload_cast<>(&ThisClass::cloneRef, py::const), "Clone the data as Ref<IData>")
    .def("cloneEmptyRef", py::overload_cast<>(&ThisClass::cloneEmptyRef, py::const), "Clone empty data")
    .def("cloneTrueRef", py::overload_cast<>(&ThisClass::cloneTrueRef, py::const), "Clone as typed ref")
    .def("cloneTrueEmptyRef", py::overload_cast<>(&ThisClass::cloneTrueEmptyRef, py::const), "Clone empty as typed ref")
    .def("storageTypeInfo", &ThisClass::storageTypeInfo, "Return storage type info")
    .def("fillDefault", &ThisClass::fillDefault, "Fill with default values")
    .def("setName", &ThisClass::setName, py::arg("name"), "Set internal name")
    .def("shape", &ThisClass::shape, "Get array shape")
    .def("setShape", &ThisClass::setShape, py::arg("new_shape"), "Set array shape")
    .def("allocationInfo", &ThisClass::allocationInfo, "Get allocation info")
    .def("setAllocationInfo", &ThisClass::setAllocationInfo, py::arg("v"), "Set allocation info")
    .def("copy", py::overload_cast<const IData*>(&ThisClass::copy),
         py::arg("data"), "Copy from another IData")
    .def("copy", py::overload_cast<const ThisClass*>(&ThisClass::copy),
         py::arg("data"), "Copy from another ArrayDataT")
    .def("swapValues", py::overload_cast<IData*>(&ThisClass::swapValues),
         py::arg("data"), "Swap values with another IData")
    .def("swapValues", py::overload_cast<ThisClass*>(&ThisClass::swapValues),
         py::arg("data"), "Swap values with another ArrayDataT")
    .def("visit", py::overload_cast<IArrayDataVisitor*>(&ThisClass::visit),
         py::arg("visitor"), "Visit with array data visitor")
    .def("visit", py::overload_cast<IDataVisitor*>(&ThisClass::visit),
         py::arg("visitor"), "Visit with generic data visitor")
    .def("__repr__", [type_name](const ThisClass& self) {
      std::ostringstream ss;
      ss << type_name << "ArrayData(size=" << self.view().size() << ")";
      return ss.str();
    }, "String representation")
    ;

  // Type alias for the typed interface
  std::string interface_name = type_name + "ArrayDataInterface";
  py::class_<DataInterfaceType, IData, std::unique_ptr<DataInterfaceType, py::deleter_base>>(m, interface_name.c_str())
    .def("view", py::overload_cast<>(&DataInterfaceType::view), "Get mutable view on data")
    .def("view", py::overload_cast<const DataInterfaceType*>(&DataInterfaceType::view, py::const), "Get const view on data")
    .def("cloneTrueRef", py::overload_cast<>(&DataInterfaceType::cloneTrueRef, py::const), "Clone as typed ref")
    .def("cloneTrueEmptyRef", py::overload_cast<>(&DataInterfaceType::cloneTrueEmptyRef, py::const), "Clone empty as typed ref")
    ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Full implementation for common types
/*---------------------------------------------------------------------------*/

void register_ArrayData(py::module_& m)
{
  // Register base types first
  register_EDataType(m);
  register_ArrayShape(m);
  register_IData(m);
  register_IArrayData(m);

  // Instantiate templates for common data types
  register_ArrayDataT_impl<Byte>(m, "Byte");
  register_ArrayDataT_impl<Int16>(m, "Int16");
  register_ArrayDataT_impl<Int32>(m, "Int32");
  register_ArrayDataT_impl<Int64>(m, "Int64");
  register_ArrayDataT_impl<Real>(m, "Real");
  register_ArrayDataT_impl<Real2>(m, "Real2");
  register_ArrayDataT_impl<Real3>(m, "Real3");
  register_ArrayDataT_impl<Real2x2>(m, "Real2x2");
  register_ArrayDataT_impl<Real3x3>(m, "Real3x3");
}

} // namespace Arcane::Pybind11

PYBIND11_MODULE(ArcaneArrayData, mod)
{
  mod.doc() = "Arcane ArrayData bindings for Python";

  py::module_ core = mod.def_submodule("core", "Arcane core bindings");

  Arcane::Pybind11::register_ArrayData(core);
}
