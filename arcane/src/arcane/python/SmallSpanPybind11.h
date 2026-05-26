// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* SmallSpanPybind11.h                                        (C) 2000-2026 */
/*                                                                           */
/* pybind11 wrapper for SmallSpan with zero-copy numpy array access.         */
/*---------------------------------------------------------------------------*/

#ifndef ARCANE_PYBIND11_SMALLSPAN_H
#define ARCANE_PYBIND11_SMALLSPAN_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "arccore/base/Span.h"

#include "arcane/utils/UtilsTypes.h"

namespace py = pybind11;

namespace Arcane::Pybind11
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Zero-copy conversion: SmallSpan<T> <-> py::array
/*---------------------------------------------------------------------------*/

namespace detail
{
  /*!\brief Convert SmallSpan<T> to py::array (zero-copy with base).
   *
   * Creates a numpy array that directly references the SmallSpan's underlying
   * data without copying. The array is read-only.
   *
   * \param span  The SmallSpan to convert.
   * \param base  Python object that owns/keeps alive the container holding the
   *              span's data. This is passed as the numpy array's base object,
   *              preventing the data from being freed while the numpy array exists.
   *
   * \returns A numpy array sharing the same data as the span (zero-copy).
   *
   * \warning The data must remain valid for as long as the returned array is used.
   */
  template<typename T>
  py::array_t<T> smallSpanToNumpyZeroCopy(const SmallSpan<T>& span, py::handle base)
  {
    std::size_t n = span.size();
    if (n == 0) {
      return py::array_t<T>();
    }

    // Use pybind11's array constructor that wraps existing memory:
    // array(void *data, dtype dt, handle base, shape, strides, offset)
    py::array_t<T> result(
        static_cast<const void*>(span.data()),  // existing data pointer
        py::dtype::of<T>(),                      // numpy dtype
        base,                                    // base keeps data alive (must not be none)
        {n}                                      // shape
    );

    // Mark as read-only
    result.mutable_flags() &= ~py::array::f_writeable;

    return result;
  }

  /*!\brief Convert SmallSpan<T> to py::array (copies data).
   *
   * Convenience overload that creates a numpy array owning its own copy.
   */
  template<typename T>
  py::array_t<T> smallSpanToNumpyCopy(const SmallSpan<T>& span)
  {
    std::size_t n = span.size();
    if (n == 0) {
      return py::array_t<T>();
    }

    py::array_t<T> result({n});
    py::buffer_info buf = result.request();
    std::copy(span.begin(), span.end(), static_cast<T*>(buf.ptr));
    return result;
  }

  /*!\brief Convert py::array_t<T> to SmallSpan<T> (zero-copy, creates view).
   *
   * Creates a SmallSpan that views the numpy array data without copying.
   * Modifications through the SmallSpan will affect the numpy array.
   *
   * \note The caller must keep the py::array_t alive while the returned
   * SmallSpan is in use.
   */
  template<typename T>
  SmallSpan<T> numpyToSmallSpanView(py::array_t<T> py_array)
  {
    py::buffer_info info = py_array.request();
    if (info.ndim != 1) {
      throw std::runtime_error("numpy array must be 1D for SmallSpan conversion");
    }
    Int32 size = static_cast<Int32>(info.shape[0]);
    return SmallSpan<T>(static_cast<T*>(info.ptr), size);
  }

  /*!\brief Convert py::buffer_info to SmallSpan<T> (zero-copy view).
   *
   * Lower-level version that takes a buffer_info directly.
   * The py::array_t must remain alive while the SmallSpan is used.
   */
  template<typename T>
  SmallSpan<T> numpyBufferToSmallSpanView(py::buffer_info info)
  {
    if (info.ndim != 1) {
      throw std::runtime_error("numpy buffer must be 1D for SmallSpan conversion");
    }
    if (info.size == 0) {
      return SmallSpan<T>();
    }
    Int32 size = static_cast<Int32>(info.shape[0]);
    return SmallSpan<T>(static_cast<T*>(info.ptr), size);
  }
} // namespace detail

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Expose SmallSpan<T, Extent> as a Python class with numpy interoperability
/*---------------------------------------------------------------------------*/

namespace detail
{
  template<typename T, Int32 Extent>
  void register_SmallSpan(py::module_& m, const std::string& type_name)
  {
    using SpanType = SmallSpan<T, Extent>;

    std::string class_name = type_name + "SmallSpan";

    py::class_<SpanType>(m, class_name.c_str())
      .def(py::init<>(), "Default constructor (empty span)")
      .def(py::init<Int32>(), py::arg("size"), "Constructor with size (creates empty span)")
      .def(py::init<T*, Int32>(), py::arg("ptr"), py::arg("size"),
           "Constructor from raw pointer and size")
      .def("size", &SpanType::size, "Number of elements in the span")
      .def("length", &SpanType::length, "Alias for size()")
      .def("size_bytes", &SpanType::sizeBytes, "Size in bytes")
      .def("empty", &SpanType::empty, "Check if span is empty")
      .def("data", [](SpanType& self) -> T* { return self.data(); },
           py::return_value_policy::reference,
           "Raw pointer to data (no ownership transfer)")
      .def("__len__", &SpanType::size, "Number of elements")
      .def("__getitem__", [](SpanType& self, std::size_t idx) -> T {
        return self[static_cast<Int32>(idx)];
      }, py::arg("idx"), "Get element at index")
      .def("__setitem__", [](SpanType& self, std::size_t idx, T val) {
        self[static_cast<Int32>(idx)] = val;
      }, py::arg("idx"), py::arg("val"), "Set element at index")
      .def("__iter__", [](SpanType& self) {
        return py::iterator(self.begin());
      }, "Iterate over elements")
      .def("subspan", py::overload_cast<Int32, Int32>(&SpanType::subspan, py::const),
           py::arg("begin"), py::arg("size"),
           "Create a subspan")
      .def("sub_span", py::overload_cast<Int32, Int32>(&SpanType::subSpan, py::const),
           py::arg("begin"), py::arg("size"),
           "Create a sub-span (alias)")
      .def("sub_part", py::overload_cast<Int32, Int32>(&SpanType::subPart, py::const),
           py::arg("begin"), py::arg("size"),
           "Create a sub-part (alias)")
      .def("fill", py::overload_cast<T>(&SpanType::fill),
           py::arg("value"), "Fill all elements with value")
      .def("__repr__", [type_name](const SpanType& self) {
        std::ostringstream ss;
        ss << type_name << "SmallSpan(size=" << self.size() << ")";
        return ss.str();
      }, "String representation")
    ;
  }
} // namespace detail

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Registration helper
/*---------------------------------------------------------------------------*/

inline void register_SmallSpanTypes(py::module_& m)
{
  detail::register_SmallSpan<std::byte>(m, "Byte");
  detail::register_SmallSpan<UChar>(m, "UChar");
  detail::register_SmallSpan<Int16>(m, "Int16");
  detail::register_SmallSpan<Int32>(m, "Int32");
  detail::register_SmallSpan<Int64>(m, "Int64");
  detail::register_SmallSpan<Integer>(m, "Integer");
  detail::register_SmallSpan<Real>(m, "Real");
  detail::register_SmallSpan<bool>(m, "Bool");
  detail::register_SmallSpan<Real2>(m, "Real2");
  detail::register_SmallSpan<Real3>(m, "Real3");
  detail::register_SmallSpan<Real2x2>(m, "Real2x2");
  detail::register_SmallSpan<Real3x3>(m, "Real3x3");
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
// Explicit converter functions for numpy <-> Arcane array views
/*---------------------------------------------------------------------------*/

namespace detail
{
  /*!\brief Convert py::array_t<T> to ArrayView<T> (zero-copy).
   *
   * Creates an ArrayView that views the numpy array data without copying.
   * The py::array_t must remain alive while the ArrayView is used.
   */
  template<typename T>
  ArrayView<T> numpyToArrayView(py::array_t<T> py_array)
  {
    py::buffer_info info = py_array.request();
    if (info.ndim != 1) {
      throw std::runtime_error("numpy array must be 1D for ArrayView conversion");
    }
    Int32 size = static_cast<Int32>(info.shape[0]);
    return ArrayView<T>(static_cast<T*>(info.ptr), size);
  }

  /*!\brief Convert py::buffer_info to ArrayView<T> (zero-copy).
   */
  template<typename T>
  ArrayView<T> numpyBufferToArrayView(py::buffer_info info)
  {
    if (info.ndim != 1) {
      throw std::runtime_error("numpy buffer must be 1D for ArrayView conversion");
    }
    if (info.size == 0) {
      return ArrayView<T>();
    }
    Int32 size = static_cast<Int32>(info.shape[0]);
    return ArrayView<T>(static_cast<T*>(info.ptr), size);
  }

  /*!\brief Convert ArrayView<T> to py::array_t<T> (copies data).
   */
  template<typename T>
  py::array_t<T> arrayViewToNumpy(ArrayView<T> view)
  {
    py::array_t<T> result(view.size());
    py::buffer_info buf = result.request();
    T* ptr = static_cast<T*>(buf.ptr);
    std::copy(view.begin(), view.end(), ptr);
    return result;
  }
} // namespace detail

} // namespace Arcane::Pybind11

#endif
