#ifndef LIB_ARCANE_ALINA_H
#define LIB_ARCANE_ALINA_H

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   lib/amgcl.h
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  C wrapper interface to amgcl.
 */

#ifdef WIN32
#define STDCALL __cdecl
#define ARCANE_ALINA_LIB_EXPORT_MACRO __declspec(dllexport)
#define ARCANE_ALINA_LIB_IMPORT_MACRO __declspec(dllimport)
#else
#  define STDCALL
#define ARCANE_ALINA_LIB_EXPORT_MACRO __attribute__ ((visibility("default")))
#define ARCANE_ALINA_LIB_IMPORT_MACRO __attribute__ ((visibility("default")))
#endif

#ifdef ARCANE_COMPONENT_arcane_alina_lib
#define ARCANE_ALINA_LIB_EXPORT ARCANE_ALINA_LIB_EXPORT_MACRO
#else
#define ARCANE_ALINA_LIB_EXPORT ARCANE_ALINA_LIB_IMPORT_MACRO
#endif

#include "arcane/alina/AlinaGlobal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void* amgclHandle;

// Create parameter list.
amgclHandle ARCANE_ALINA_LIB_EXPORT STDCALL ARCANE_ALINA_params_create();

// Set integer parameter in a parameter list.
void ARCANE_ALINA_EXPORT STDCALL ARCANE_ALINA_params_seti(amgclHandle prm, const char* name, int value);

// Set floating point parameter in a parameter list.
void ARCANE_ALINA_EXPORT STDCALL ARCANE_ALINA_params_setf(amgclHandle prm, const char* name, float value);

// Set floating point parameter in a parameter list.
void ARCANE_ALINA_EXPORT STDCALL ARCANE_ALINA_params_sets(amgclHandle prm, const char* name, const char* value);

// Read parameters from a JSON file
void ARCANE_ALINA_EXPORT STDCALL ARCANE_ALINA_params_read_json(amgclHandle prm, const char* fname);

// Destroy parameter list.
void ARCANE_ALINA_EXPORT STDCALL ARCANE_ALINA_params_destroy(amgclHandle prm);

// Create AMG preconditioner.
amgclHandle ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_precond_create(int n,
                     const int* ptr,
                     const int* col,
                     const double* val,
                     amgclHandle parameters);

// Apply AMG preconditioner (x = M^(-1) * rhs).
void ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_precond_apply(amgclHandle amg, const double* rhs, double* x);

// Printout preconditioner structure
void ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_precond_report(amgclHandle amg);

// Destroy AMG preconditioner
void ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_precond_destroy(amgclHandle amg);

// Create iterative solver preconditioned by AMG.
amgclHandle ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_solver_create(int n,
                    const int* ptr,
                    const int* col,
                    const double* val,
                    amgclHandle parameters);

// Convergence info
struct ARCANE_ALINA_LIB_EXPORT conv_info
{
  int iterations;
  double residual;
};

// Solve the problem for the given right-hand side.
conv_info ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_solver_solve(amgclHandle solver,
                   double const* rhs,
                   double* x);

// Solve the problem for the given matrix and the right-hand side.
conv_info ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_solver_solve_mtx(amgclHandle solver,
                       int const* A_ptr,
                       int const* A_col,
                       double const* A_val,
                       double const* rhs,
                       double* x);

// Printout solver structure
void ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_solver_report(amgclHandle solver);

// Destroy iterative solver.
void ARCANE_ALINA_EXPORT STDCALL
ARCANE_ALINA_solver_destroy(amgclHandle solver);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
