// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* ArrayDataPybind11.h                                       (C) 2000-2026 */
/*                                                                           */
/* pybind11 wrapper for ArrayDataT and related classes.                      */
/*---------------------------------------------------------------------------*/

#ifndef ARCANE_PYBIND11_ARRAYDATA_H
#define ARCANE_PYBIND11_ARRAYDATA_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "arcane/core/IData.h"
#include "arcane/core/IDataVisitor.h"
#include "arcane/core/ISerializer.h"
#include "arcane/core/IDataOperation.h"
#include "arcane/core/IHashAlgorithm.h"
#include "arcane/core/datatype/DataAllocationInfo.h"
#include "arcane/utils/Array.h"
#include "arcane/utils/ArrayShape.h"
#include "arcane/utils/Ref.h"
#include "arcane/impl/internal/ArrayData.h"

#include "ArcanePybind11Global.h"

namespace py = pybind11;

namespace Arcane::Pybind11
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Wrapper pour eDataType
/*---------------------------------------------------------------------------*/

inline void register_EDataType(py::module_& m)
{
  py::enum_<eDataType>(m, "EDataType")
    .value("DT_None", eDataType::DT_None)
    .value("DT_Byte", eDataType::DT_Byte)
    .value("DT_Int16", eDataType::DT_Int16)
    .value("DT_Int32", eDataType::DT_Int32)
    .value("DT_Int64", eDataType::DT_Int64)
    .value("DT_Real", eDataType::DT_Real)
    .value("DT_Real2", eDataType::DT_Real2)
    .value("DT_Real3", eDataType::DT_Real3)
    .value("DT_Real2x2", eDataType::DT_Real2x2)
    .value("DT_Real3x3", eDataType::DT_Real3x3)
    .export_values();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Wrapper pour ArrayShape
/*---------------------------------------------------------------------------*/

inline void register_ArrayShape(py::module_& m)
{
  py::class_<ArrayShape>(m, "ArrayShape")
    .def(py::init<>(), "Default constructor")
    .def(py::init<const Int64Array&>(), "Constructor from Int64Array")
    .def("size", &ArrayShape::size, "Return number of dimensions")
    .def("totalNbElement", &ArrayShape::totalNbElement, "Return total number of elements")
    .def("__len__", &ArrayShape::size, "Number of dimensions")
    .def("__getitem__", [](const ArrayShape& s, std::size_t idx) -> Int64 {
      return s[idx];
    }, py::arg("idx"), "Get dimension size at index")
    .def("__setitem__", [](ArrayShape& s, std::size_t idx, Int64 val) {
      s[idx] = val;
    }, py::arg("idx"), py::arg("val"), "Set dimension size at index")
    .def("__iter__", [](ArrayShape& s) -> py::iterator {
      return py::iterator(s.begin());
    }, "Iterate over dimensions")
    .def("__repr__", [](const ArrayShape& s) {
      std::ostringstream ss;
      ss << "ArrayShape(";
      for (std::size_t i = 0; i < s.size(); ++i) {
        if (i > 0) ss << ", ";
        ss << s[i];
      }
      ss << ")";
      return ss.str();
    }, "String representation")
    .def("__eq__", &ArrayShape::operator==, py::arg("other"), "Equality comparison")
    ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Wrapper pour IData (interface de base)
/*---------------------------------------------------------------------------*/

inline void register_IData(py::module_& m)
{
  py::class_<IData, std::unique_ptr<IData, py::deleter_base>>(m, "IData")
    .def("dimension", &IData::dimension, "Return data dimension (0=scalar, 1=array)")
    .def("multiTag", &IData::multiTag, "Return multi-tag value")
    .def("dataType", &IData::dataType, "Return data type")
    .def("cloneRef", &IData::cloneRef, "Clone the data")
    .def("cloneEmptyRef", &IData::cloneEmptyRef, "Clone empty data")
    .def("storageTypeInfo", &IData::storageTypeInfo, "Return storage type info")
    .def("resize", &IData::resize, py::arg("new_size"), "Resize 1D array data")
    .def("fillDefault", &IData::fillDefault, "Fill with default values")
    .def("setName", &IData::setName, py::arg("name"), "Set internal name")
    .def("shape", &IData::shape, "Return array shape")
    .def("setShape", &IData::setShape, py::arg("new_shape"), "Set array shape")
    .def("allocationInfo", &IData::allocationInfo, "Return allocation info")
    .def("setAllocationInfo", &IData::setAllocationInfo, py::arg("v"), "Set allocation info")
    .def("copy", &IData::copy, py::arg("data"), "Copy data from another IData")
    .def("swapValues", &IData::swapValues, py::arg("data"), "Swap values with another IData")
    ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Wrapper pour IArrayData
/*---------------------------------------------------------------------------*/

inline void register_IArrayData(py::module_& m)
{
  py::class_<IArrayData, IData, std::unique_ptr<IArrayData, py::deleter_base>>(m, "IArrayData")
    ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Template instantiation pour ArrayDataT<DataType>
/*---------------------------------------------------------------------------*/

template<typename DataType>
void register_ArrayDataT(py::module_& m, const std::string& type_name)
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
    .def("serialize", py::overload_cast<ISerializer*, IDataOperation*>(&ThisClass::serialize),
         py::arg("sbuf"), py::arg("operation"), "Serialize data")
    .def("serialize", py::overload_cast<ISerializer*, Int32ConstArrayView, IDataOperation*>(&ThisClass::serialize),
         py::arg("sbuf"), py::arg("ids"), py::arg("operation"), "Serialize subset of data")
    .def("fillDefault", &ThisClass::fillDefault, "Fill with default values")
    .def("setName", &ThisClass::setName, py::arg("name"), "Set internal name")
    .def("createSerializedDataRef", py::overload_cast<bool>(&ThisClass::createSerializedDataRef, py::const),
         py::arg("use_basic_type"), "Create serialized data")
    .def("allocateBufferForSerializedData", &ThisClass::allocateBufferForSerializedData,
         py::arg("sdata"), "Allocate buffer for serialized data")
    .def("assignSerializedData", &ThisClass::assignSerializedData,
         py::arg("sdata"), "Assign serialized data")
    .def("copy", py::overload_cast<const IData*>(&ThisClass::copy),
         py::arg("data"), "Copy from another IData")
    .def("copy", py::overload_cast<const ThisClass*>(&ThisClass::copy),
         py::arg("data"), "Copy from another ArrayDataT")
    .def("swapValues", py::overload_cast<IData*>(&ThisClass::swapValues),
         py::arg("data"), "Swap values with another IData")
    .def("swapValues", py::overload_cast<ThisClass*>(&ThisClass::swapValues),
         py::arg("data"), "Swap values with another ArrayDataT")
    .def("computeHash", &ThisClass::computeHash,
         py::arg("algo"), py::arg("output"), "Compute hash")
    .def("shape", &ThisClass::shape, "Get array shape")
    .def("setShape", &ThisClass::setShape, py::arg("new_shape"), "Set array shape")
    .def("allocationInfo", &ThisClass::allocationInfo, "Get allocation info")
    .def("setAllocationInfo", &ThisClass::setAllocationInfo, py::arg("v"), "Set allocation info")
    .def("visit", py::overload_cast<IArrayDataVisitor*>(&ThisClass::visit),
         py::arg("visitor"), "Visit with array data visitor")
    .def("visit", py::overload_cast<IDataVisitor*>(&ThisClass::visit),
         py::arg("visitor"), "Visit with generic data visitor")
    .def("_internal", &ThisClass::_internal, py::return_value_policy::reference,
         "Get internal interface (raw pointer, caller must not delete)")
    .def("_commonInternal", &ThisClass::_commonInternal, py::return_value_policy::reference,
         "Get common internal interface (raw pointer, caller must not delete)")
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
    .def("_internal", &DataInterfaceType::_internal, py::return_value_policy::reference,
         "Get internal interface")
    ;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Registration function for common data types
/*---------------------------------------------------------------------------*/

inline void register_ArrayData(py::module_& m)
{
  // Register base types first
  register_EDataType(m);
  register_ArrayShape(m);
  register_IData(m);
  register_IArrayData(m);

  // Instantiate templates for common data types
  register_ArrayDataT<Byte>(m, "Byte");
  register_ArrayData<Int16>(m, "Int16");
  register_ArrayData<Int32>(m, "Int32");
  register_ArrayData<Int64>(m, "Int64");
  register_ArrayData<Real>(m, "Real");
  register_ArrayData<Real2>(m, "Real2");
  register_ArrayData<Real3>(m, "Real3");
  register_ArrayData<Real2x2>(m, "Real2x2");
  register_ArrayData<Real3x3>(m, "Real3x3");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Helper function to convert Arcane ArrayView to/from numpy (optional)
/*---------------------------------------------------------------------------*/

namespace detail
{
  template<typename T>
  py::array_t<T> arrayViewToNumpy(ArrayView<T> view)
  {
    py::array_t<T> py_array(view.size());
    py::buffer_info buf = py_array.request();
    T* ptr = static_cast<T*>(buf.ptr);
    std::copy(view.begin(), view.end(), ptr);
    return py_array;
  }

  template<typename T>
  ArrayView<T> numpyToArrayView(py::array_t<T> py_array)
  {
    py::buffer_info info = py_array.request();
    if (info.ndim != 1)
      ARCANE_THROW(Argexption, "numpy array must be 1D");
    return ArrayView<T>(static_cast<T*>(info.ptr), info.size);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Pybind11

#endif
