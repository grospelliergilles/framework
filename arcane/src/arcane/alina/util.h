// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Util.h                                                      (C) 2026-2026 */
/*                                                                           */
/* Various utilities.                                                        */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_UTIL_H
#define ARCANE_ALINA_UTIL_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*
 * This file is based on the work on AMGCL library (version march 2026)
 * which can be found at https://github.com/ddemidov/amgcl.
 *
 * Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>
 * SPDX-License-Identifier: MIT
 */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#pragma GCC diagnostic ignored "-Wconversion"

#include "arcane/alina/AlinaGlobal.h"

#include "arccore/base/Ref.h"
#include "arccore/base/FatalErrorException.h"

#include <set>
#include <complex>
#include <cstddef>
#include <tuple>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
//! Result of a solving.
class ARCANE_ALINA_EXPORT SolverResult
{
 public:

  SolverResult() = default;
  SolverResult(const std::tuple<size_t, double>& v)
  : m_nb_iteration(get<0>(v))
  , m_residual(get<1>(v))
  {}
  SolverResult(const std::tuple<size_t, float>& v)
  : m_nb_iteration(get<0>(v))
  , m_residual(get<1>(v))
  {}
  SolverResult(size_t nb_iteration, double residual)
  : m_nb_iteration(nb_iteration)
  , m_residual(residual)
  {}

  operator std::tuple<size_t, double>() const { return { m_nb_iteration, m_residual }; }

 public:

  constexpr Int32 nbIteration() const { return static_cast<Int32>(m_nb_iteration); }
  constexpr double residual() const { return m_residual; }

 private:

  size_t m_nb_iteration = 0;
  double m_residual = 0.0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// If asked explicitly, or if boost is available, enable
// using boost::propert_tree::ptree as parameters:
#include <boost/property_tree/ptree.hpp>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Class to wrap 'boost::property_tree::ptree' to ease removing it
class ARCANE_ALINA_EXPORT PropertyTree
{
 public:

  using BoostPTree = boost::property_tree::ptree;

 public:

  PropertyTree();
  PropertyTree(const PropertyTree& rhs);
  explicit PropertyTree(const BoostPTree& x);
  ~PropertyTree();

 public:

  template <typename DataType> auto
  get(const char* param_type, const DataType& default_value) const
  {
    return m_property_tree->get(param_type, default_value);
  }
  template <typename DataType> void
  put(const std::string& path, const DataType& value)
  {
    m_property_tree->put(path, value);
  }
  PropertyTree get_child_empty(const std::string& path) const;
  bool erase(const char* name);
  size_t count(const char* name) const;

 public:

  void read_json(const std::string& filename);

 public:

  //TODO: need remove
  const BoostPTree& toBoostPTree() const { return *m_property_tree; }
  //TODO: need remove
  BoostPTree& toBoostPTree() { return *m_property_tree; }

 private:

  BoostPTree* m_property_tree = nullptr;
  bool m_is_own = false;
};

namespace detail
{

  struct empty_params
  {
    empty_params() {}

    empty_params(const PropertyTree& ap)
    {
      const boost::property_tree::ptree& p = ap.toBoostPTree();
      for (const auto& v : p) {
        std::cerr << "Alina: unknown parameter " << v.first << "\n";
      }
    }
    void get(boost::property_tree::ptree&, const std::string&) const {}
  };

} // namespace detail

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

namespace boost::property_tree::json_parser
{

inline void
read_json(const std::string& filename, Arcane::Alina::PropertyTree& prm)
{
  prm.read_json(filename);
}

} // namespace boost::property_tree::json_parser

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/

#include <arcane/alina/ScopedStreamModifier.h>

/*!
 * \brief Performance measurement macros.
 *
 * If ARCANE_ALINA_PROFILING macro is defined at compilation, then ARCANE_ALINA_TIC(name) and
 * ARCANE_ALINA_TOC(name) macros correspond to prof.tic(name) and prof.toc(name).
 * Arcane::Alina::prof should be an instance of Arcane::Alina::profiler defined in a user
 * code similar to:
 * \code
 * namespace Arcane::Alina { profiler prof; }
 * \endcode
 * If ARCANE_ALINA_PROFILING is undefined, then ARCANE_ALINA_TIC and ARCANE_ALINA_TOC are noop macros.
 */
#ifdef ARCANE_ALINA_PROFILING
#  if !defined(ARCANE_ALINA_TIC) || !defined(ARCANE_ALINA_TOC)
#    include <arcane/alina/Profiler.h>
#    define ARCANE_ALINA_TIC(name):: Arcane::Alina::Profiler::globalTic(name);
#    define ARCANE_ALINA_TOC(name) ::Arcane::Alina::Profiler::globalToc(name);
#  endif
#endif

#ifndef ARCANE_ALINA_TIC
#  define ARCANE_ALINA_TIC(name)
#endif
#ifndef ARCANE_ALINA_TOC
#  define ARCANE_ALINA_TOC(name)
#endif

#define ARCANE_ALINA_DEBUG_SHOW(x)                                                    \
    std::cout << std::setw(20) << #x << ": "                                   \
              << std::setw(15) << std::setprecision(8) << std::scientific      \
              << (x) << std::endl

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Throws \p message if \p condition is not true.
template <class Condition, class Message>
void precondition(const Condition& condition, const Message& message)
{
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4800)
#endif
  if (!condition)
    throw std::runtime_error(message);
#ifdef _MSC_VER
#pragma warning(pop)
#endif
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#define ARCANE_ALINA_PARAMS_IMPORT_VALUE(p, name) \
  name(p.get(#name, params().name))

#define ARCANE_ALINA_PARAMS_IMPORT_CHILD(p, name) \
  name(p.get_child_empty(#name))

#define ARCANE_ALINA_PARAMS_EXPORT_VALUE(p, path, name) \
  p.put(std::string(path) + #name, name)

namespace detail
{

  template <typename T>
  inline void params_export_child(PropertyTree& ap,
                                  const std::string& path,
                                  const char* name, const T& obj)
  {
    boost::property_tree::ptree& p = ap.toBoostPTree();
    obj.get(p, std::string(path) + name + ".");
  }

  template <>
  inline void params_export_child(PropertyTree& ap,
                                  const std::string& path, const char* name,
                                  const boost::property_tree::ptree& obj)
  {
    boost::property_tree::ptree& p = ap.toBoostPTree();
    p.add_child(std::string(path) + name, obj);
  }

} // namespace detail

#define ARCANE_ALINA_PARAMS_EXPORT_CHILD(p, path, name) \
  ::Arcane::Alina::detail::params_export_child(p, path, #name, name)

// Missing parameter action
#ifndef ARCANE_ALINA_PARAM_MISSING
#define ARCANE_ALINA_PARAM_MISSING(name) (void)0
#endif

// Unknown parameter action
#ifndef ARCANE_ALINA_PARAM_UNKNOWN
#define ARCANE_ALINA_PARAM_UNKNOWN(name) \
  std::cerr << "AMGCL WARNING: unknown parameter " << name << std::endl
#endif

extern "C++" ARCANE_ALINA_EXPORT void
check_params(const PropertyTree& ptree, const std::set<std::string>& names);

extern "C++" ARCANE_ALINA_EXPORT void
check_params(const PropertyTree& ptree,
             const std::set<std::string>& names,
             const std::set<std::string>& opt_names);

// Put parameter in form "key=value" into a boost::property_tree::ptree
extern "C++" ARCANE_ALINA_EXPORT void
put(PropertyTree& ptree, const std::string& param);

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Iterator range
template <class Iterator>
class iterator_range
{
 public:

  typedef Iterator iterator;
  typedef Iterator const_iterator;
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
  typedef typename std::iterator_traits<Iterator>::reference reference;

  iterator_range(Iterator b, Iterator e)
  : b(b)
  , e(e)
  {}

  ptrdiff_t size() const
  {
    return std::distance(b, e);
  }

  Iterator begin() const
  {
    return b;
  }

  Iterator end() const
  {
    return e;
  }

  reference operator[](size_t i) const
  {
    return b[i];
  }

 private:

  Iterator b, e;
};

template <class Iterator>
iterator_range<Iterator> make_iterator_range(Iterator b, Iterator e)
{
  return iterator_range<Iterator>(b, e);
}

// N-dimensional dense matrix
template <class T, int N>
class multi_array
{
  static_assert(N > 0, "Wrong number of dimensions");

 public:

  template <class... I>
  multi_array(I... n)
  {
    static_assert(sizeof...(I) == N, "Wrong number of dimensions");
    buf.resize(init(n...));
  }

  size_t size() const
  {
    return buf.size();
  }

  int stride(int i) const
  {
    return strides[i];
  }

  template <class... I>
  T operator()(I... i) const
  {
    static_assert(sizeof...(I) == N, "Wrong number of indices");
    return buf[index(i...)];
  }

  template <class... I>
  T& operator()(I... i)
  {
    static_assert(sizeof...(I) == N, "Wrong number of indices");
    return buf[index(i...)];
  }

  const T* data() const
  {
    return buf.data();
  }

  T* data()
  {
    return buf.data();
  }

 private:

  std::array<int, N> strides;
  std::vector<T> buf;

  template <class... I>
  int index(int i, I... tail) const
  {
    return strides[N - sizeof...(I) - 1] * i + index(tail...);
  }

  int index(int i) const
  {
    return strides[N - 1] * i;
  }

  template <class... I>
  int init(int i, I... tail)
  {
    int size = init(tail...);
    strides[N - sizeof...(I) - 1] = size;
    return i * size;
  }

  int init(int i)
  {
    strides[N - 1] = 1;
    return i;
  }
};

template <class T>
class circular_buffer
{
 public:

  circular_buffer(size_t n)
  : start(0)
  {
    buf.reserve(n);
  }

  size_t size() const
  {
    return buf.size();
  }

  void push_back(const T& v)
  {
    if (buf.size() < buf.capacity()) {
      buf.push_back(v);
    }
    else {
      buf[start] = v;
      start = (start + 1) % buf.capacity();
    }
  }

  const T& operator[](size_t i) const
  {
    return buf[(start + i) % buf.capacity()];
  }

  T& operator[](size_t i)
  {
    return buf[(start + i) % buf.capacity()];
  }

  void clear()
  {
    buf.clear();
    start = 0;
  }

 private:

  size_t start;
  std::vector<T> buf;
};

namespace detail
{

  template <class T>
  T eps(size_t n)
  {
    return 2 * std::numeric_limits<T>::epsilon() * n;
  }

} // namespace detail

template <class T> struct is_complex : std::false_type
{};
template <class T> struct is_complex<std::complex<T>> : std::true_type
{};

inline std::string human_readable_memory(size_t bytes)
{
  static const char* suffix[] = { "B", "K", "M", "G", "T" };

  int i = 0;
  double m = static_cast<double>(bytes);
  for (; i < 4 && m >= 1024.0; ++i, m /= 1024.0)
    ;

  std::ostringstream s;
  s << std::fixed << std::setprecision(2) << m << " " << suffix[i];
  return s.str();
}

namespace detail
{

  class non_copyable
  {
   protected:

    non_copyable() = default;
    ~non_copyable() = default;

    non_copyable(non_copyable const&) = delete;
    void operator=(non_copyable const& x) = delete;
  };

  /*!
   * \brief  Sort row of CRS matrix by columns.
   */
  template <typename Col, typename Val>
  void sort_row(Col* col, Val* val, int n)
  {
    for (int j = 1; j < n; ++j) {
      Col c = col[j];
      Val v = val[j];

      int i = j - 1;

      while (i >= 0 && col[i] > c) {
        col[i + 1] = col[i];
        val[i + 1] = val[i];
        i--;
      }

      col[i + 1] = c;
      val[i + 1] = v;
    }
  }

} // namespace detail

namespace error
{

  struct empty_level
  {};

} // namespace error
} // namespace Arcane::Alina

namespace std
{

// Read pointers from input streams.
// This allows to exchange pointers through boost::property_tree::ptree.
template <class T>
inline istream& operator>>(istream& is, T*& ptr)
{
  Arcane::Alina::ScopedStreamModifier ss(is);

  size_t val;
  is >> std::hex >> val;

  ptr = reinterpret_cast<T*>(val);

  return is;
}

} // namespace std

#endif
