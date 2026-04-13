// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* AlinaUtils.cc                                               (C) 2000-2026 */
/*                                                                           */
/* Classes utilitaires.                                                      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/ArcaneGlobal.h"

#include "arcane/alina/util.h"

#include <boost/property_tree/json_parser.hpp>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace detail
{
  inline const boost::property_tree::ptree& empty_ptree()
  {
    static const boost::property_tree::ptree p;
    return p;
  }

} // namespace detail

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PropertyTree::
PropertyTree()
: m_property_tree(new BoostPTree())
, m_is_own(true)
{
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PropertyTree::
PropertyTree(const PropertyTree& rhs)
{
  if (rhs.m_is_own) {
    m_property_tree = new BoostPTree(*rhs.m_property_tree);
    m_is_own = true;
  }
  else {
    m_property_tree = rhs.m_property_tree;
      m_is_own = false;
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PropertyTree::
PropertyTree(const BoostPTree& x)
: m_property_tree(new BoostPTree(x))
, m_is_own(true)
{}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PropertyTree::
~PropertyTree()
{
  if (m_is_own)
    delete m_property_tree;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

PropertyTree PropertyTree::
get_child_empty(const std::string& path) const
{
  const BoostPTree& child = m_property_tree->get_child(path, detail::empty_ptree());
  PropertyTree p;
  p.m_property_tree = const_cast<BoostPTree*>(&child);
  p.m_is_own = false;
  return p;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

bool PropertyTree::
erase(const char* name)
{
  return m_property_tree->erase(name);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

size_t PropertyTree::
count(const char* name) const
{
  return m_property_tree->count(name);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void PropertyTree::
read_json(const std::string& filename)
{
  boost::property_tree::ptree& p = *m_property_tree;
  boost::property_tree::json_parser::read_json(filename, p);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void check_params(const PropertyTree& ptree,
                  const std::set<std::string>& names)
{
  const boost::property_tree::ptree& p = ptree.toBoostPTree();

  for (const auto& n : names) {
    if (!p.count(n)) {
      ARCANE_ALINA_PARAM_MISSING(n);
    }
  }
  for (const auto& v : p) {
    if (!names.count(v.first)) {
      ARCANE_ALINA_PARAM_UNKNOWN(v.first);
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void check_params(const PropertyTree& ptree,
                  const std::set<std::string>& names,
                  const std::set<std::string>& opt_names)
{
  const boost::property_tree::ptree& p = ptree.toBoostPTree();

  for (const auto& n : names) {
    if (!p.count(n)) {
      ARCANE_ALINA_PARAM_MISSING(n);
    }
  }
  for (const auto& n : opt_names) {
    if (!p.count(n)) {
      ARCANE_ALINA_PARAM_MISSING(n);
    }
  }
  for (const auto& v : p) {
    if (!names.count(v.first) && !opt_names.count(v.first)) {
      ARCANE_ALINA_PARAM_UNKNOWN(v.first);
    }
  }
}

void put(PropertyTree& ptree, const std::string& param)
{
  boost::property_tree::ptree& p = ptree.toBoostPTree();
  size_t eq_pos = param.find('=');
  if (eq_pos == std::string::npos)
    throw std::invalid_argument("param in put() should have \"key=value\" format!");
  p.put(param.substr(0, eq_pos), param.substr(eq_pos + 1));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
