// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2026-2026 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* PreconditionerRuntime.h                                     (C) 2026-2026 */
/*                                                                           */
/* Runtime-configurable preconditioners.                                     */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_ALINA_PRECONDITIONERRUNTIME_H
#define ARCANE_ALINA_PRECONDITIONERRUNTIME_H
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

#include <arcane/alina/util.h>
#include <arcane/alina/solver_runtime.h>
#include <arcane/alina/CoarseningRuntime.h>
#include <arcane/alina/RelaxationRuntime.h>
#include <arcane/alina/relaxation.h>
#include <arcane/alina/DummyPreconditioner.h>
#include <arcane/alina/make_solver.h>
#include <arcane/alina/AMG.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina::runtime
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/// Preconditioner kinds.
namespace precond_class
{
  enum type
  {
    amg, ///< AMG
    relaxation, ///< Single-level relaxation
    dummy, ///< Identity matrix as preconditioner.
    nested ///< Nested solver as preconditioner.
  };

  inline std::ostream& operator<<(std::ostream& os, type p)
  {
    switch (p) {
    case amg:
      return os << "amg";
    case relaxation:
      return os << "relaxation";
    case dummy:
      return os << "dummy";
    case nested:
      return os << "nested";
    default:
      return os << "???";
    }
  }

  inline std::istream& operator>>(std::istream& in, type& p)
  {
    std::string val;
    in >> val;

    if (val == "amg")
      p = amg;
    else if (val == "relaxation")
      p = relaxation;
    else if (val == "dummy")
      p = dummy;
    else if (val == "nested")
      p = nested;
    else
      throw std::invalid_argument("Invalid preconditioner class. Valid choices are: "
                                  "amg, relaxation, dummy, nested");

    return in;
  }
} // namespace precond_class

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Runtime-configurable preconditioners.
 */
template <class Backend>
class PreconditionerRuntime
{
 public:

  typedef Backend backend_type;

  typedef typename Backend::value_type value_type;
  typedef typename Backend::matrix matrix;
  typedef typename Backend::vector vector;
  typedef typename Backend::params backend_params;

  typedef Alina::PropertyTree params;

  template <class Matrix>
  PreconditionerRuntime(const Matrix& A,
                        params prm = params(),
                        const backend_params& bprm = backend_params())
  : _class(prm.get("class", runtime::precond_class::amg))
  , handle(0)
  {
    if (!prm.erase("class"))
      ARCANE_ALINA_PARAM_MISSING("class");
    std::cout << "PreconditionerClass=" << _class << "\n";
    switch (_class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      handle = static_cast<void*>(new Precond(A, prm, bprm));
    } break;
    case precond_class::relaxation: {
      typedef Alina::relaxation::as_preconditioner<Backend, runtime::relaxation::RuntimeRelaxation>
      Precond;

      handle = static_cast<void*>(new Precond(A, prm, bprm));
    } break;
    case precond_class::dummy: {
      typedef Alina::preconditioner::DummyPreconditioner<Backend>
      Precond;

      handle = static_cast<void*>(new Precond(A, prm, bprm));
    } break;
    case precond_class::nested: {
      typedef make_solver<
      PreconditionerRuntime,
      runtime::solver::wrapper<Backend>>
      Precond;

      handle = static_cast<void*>(new Precond(A, prm, bprm));
    } break;
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

  ~PreconditionerRuntime()
  {
    switch (_class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      delete static_cast<Precond*>(handle);
    } break;
    case precond_class::relaxation: {
      typedef Alina::relaxation::as_preconditioner<Backend, runtime::relaxation::RuntimeRelaxation>
      Precond;

      delete static_cast<Precond*>(handle);
    } break;
    case precond_class::dummy: {
      typedef Alina::preconditioner::DummyPreconditioner<Backend>
      Precond;

      delete static_cast<Precond*>(handle);
    } break;
    case precond_class::nested: {
      typedef make_solver<
      PreconditionerRuntime,
      runtime::solver::wrapper<Backend>>
      Precond;

      delete static_cast<Precond*>(handle);
    } break;
    default:
      break;
    }
  }

  template <class Matrix>
  void rebuild(const Matrix& A,
               const backend_params& bprm = backend_params())
  {
    switch (_class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      static_cast<Precond*>(handle)->rebuild(A, bprm);
    } break;
    default:
      std::cerr << "rebuild is a noop unless the preconditioner is AMG" << std::endl;
      return;
    }
  }

  template <class Vec1, class Vec2>
  void apply(const Vec1& rhs, Vec2& x) const
  {
    std::cout << "ApplyPrecond class=" << _class << "\n";
    switch (_class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      static_cast<Precond*>(handle)->apply(rhs, x);
    } break;
    case precond_class::relaxation: {
      typedef Alina::relaxation::as_preconditioner<Backend, runtime::relaxation::RuntimeRelaxation>
      Precond;

      static_cast<Precond*>(handle)->apply(rhs, x);
    } break;
    case precond_class::dummy: {
      typedef Alina::preconditioner::DummyPreconditioner<Backend>
      Precond;

      static_cast<Precond*>(handle)->apply(rhs, x);
    } break;
    case precond_class::nested: {
      typedef make_solver<
      PreconditionerRuntime,
      runtime::solver::wrapper<Backend>>
      Precond;

      static_cast<Precond*>(handle)->apply(rhs, x);
    } break;
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

  std::shared_ptr<matrix> system_matrix_ptr() const
  {
    switch (_class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      return static_cast<Precond*>(handle)->system_matrix_ptr();
    }
    case precond_class::relaxation: {
      typedef Alina::relaxation::as_preconditioner<Backend, runtime::relaxation::RuntimeRelaxation>
      Precond;

      return static_cast<Precond*>(handle)->system_matrix_ptr();
    }
    case precond_class::dummy: {
      typedef Alina::preconditioner::DummyPreconditioner<Backend>
      Precond;

      return static_cast<Precond*>(handle)->system_matrix_ptr();
    }
    case precond_class::nested: {
      typedef make_solver<
      PreconditionerRuntime,
      runtime::solver::wrapper<Backend>>
      Precond;

      return static_cast<Precond*>(handle)->system_matrix_ptr();
    }
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

  const matrix& system_matrix() const
  {
    return *system_matrix_ptr();
  }

  size_t size() const
  {
    return backend::rows(system_matrix());
  }

  size_t bytes() const
  {
    switch (_class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      return backend::bytes(*static_cast<Precond*>(handle));
    }
    case precond_class::relaxation: {
      typedef Alina::relaxation::as_preconditioner<Backend, runtime::relaxation::RuntimeRelaxation>
      Precond;

      return backend::bytes(*static_cast<Precond*>(handle));
    }
    case precond_class::dummy: {
      typedef Alina::preconditioner::DummyPreconditioner<Backend>
      Precond;

      return backend::bytes(*static_cast<Precond*>(handle));
    }
    case precond_class::nested: {
      typedef make_solver<
      PreconditionerRuntime,
      runtime::solver::wrapper<Backend>>
      Precond;

      return backend::bytes(*static_cast<Precond*>(handle));
    }
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const PreconditionerRuntime& p)
  {
    switch (p._class) {
    case precond_class::amg: {
      typedef Alina::AMG<Backend, runtime::coarsening::CoarseningRuntime, runtime::relaxation::RuntimeRelaxation>
      Precond;

      return os << *static_cast<Precond*>(p.handle);
    }
    case precond_class::relaxation: {
      typedef Alina::relaxation::as_preconditioner<Backend, runtime::relaxation::RuntimeRelaxation>
      Precond;

      return os << *static_cast<Precond*>(p.handle);
    }
    case precond_class::dummy: {
      typedef Alina::preconditioner::DummyPreconditioner<Backend>
      Precond;

      return os << *static_cast<Precond*>(p.handle);
    }
    case precond_class::nested: {
      typedef make_solver<
      PreconditionerRuntime,
      runtime::solver::wrapper<Backend>>
      Precond;

      return os << *static_cast<Precond*>(p.handle);
    }
    default:
      throw std::invalid_argument("Unsupported preconditioner class");
    }
  }

 private:

  const runtime::precond_class::type _class;

  void* handle = nullptr;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina::runtime

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
