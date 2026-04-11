#ifndef ARCANE_ALINA_SOLVER_PRECOND_SIDE_HPP
#define ARCANE_ALINA_SOLVER_PRECOND_SIDE_HPP

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   alina/solver/precond_side.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Definitions and functions supporting left/right preconditioning.
 */

#include <iostream>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include <arcane/alina/ValueTypeInterface.h>

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arcane::Alina
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

enum class ePreconditionerSideType
{
  left,
  right
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline std::ostream& operator<<(std::ostream& os, ePreconditionerSideType p)
{
  switch (p) {
  case ePreconditionerSideType::left:
    return os << "left";
  case ePreconditionerSideType::right:
    return os << "right";
  default:
    return os << "???";
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

inline std::istream& operator>>(std::istream& in, ePreconditionerSideType& p)
{
  std::string val;
  in >> val;

  if (val == "left")
    p = ePreconditionerSideType::left;
  else if (val == "right")
    p = ePreconditionerSideType::right;
  else
    throw std::invalid_argument("Invalid preconditioning side. "
                                "Valid choices are: left, right.");

  return in;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Preconditioned matrix-vector product
template <class Precond, class Matrix, class VecF, class VecX, class VecT> void
preconditioner_spmv(ePreconditionerSideType pside, const Precond& P, const Matrix& A,
                    const VecF& F, VecX& X, VecT& T)
{
  typedef typename backend::value_type<Matrix>::type value;
  typedef typename math::scalar_of<value>::type scalar;

  static const scalar one = math::identity<scalar>();
  static const scalar zero = math::zero<scalar>();

  if (pside == ePreconditionerSideType::left) {
    backend::spmv(one, A, F, zero, T);
    P.apply(T, X);
  }
  else {
    P.apply(F, T);
    backend::spmv(one, A, T, zero, X);
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // namespace Arcane::Alina

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
