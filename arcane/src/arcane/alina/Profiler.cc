// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Profiler.cc                                                 (C) 2000-2026 */
/*                                                                           */
/* Classes utilitaires.                                                      */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/alina/Profiler.h"

#include "arcane/alina/perf_counter_clock.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

using Counter = perf_counter::clock;

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Profiler::
tic(const std::string& name)
{
  stack.back()->children[name].begin = Counter::current();
  stack.push_back(&stack.back()->children[name]);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

Profiler::delta_type Profiler::
toc(const std::string&)
{
  profile_unit* top = stack.back();
  stack.pop_back();

  value_type current = Counter::current();
  delta_type delta = current - top->begin;

  top->length += delta;
  root.length = current - root.begin;

  return delta;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Profiler::
reset()
{
  stack.clear();
  root.length = 0;
  root.children.clear();

  stack.push_back(&root);
  root.begin = Counter::current();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
void Profiler::
init()
{
  stack.reserve(128);
  stack.push_back(&root);
  root.begin = Counter::current();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void Profiler::
print(std::ostream& out)
{
  if (stack.back() != &root)
    out << "Warning! Profile is incomplete." << std::endl;
  ios_saver ss(out);
  root.print(out, name, 0, root.length, root.total_width(name, 0));
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
