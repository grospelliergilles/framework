﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
//-----------------------------------------------------------------------------
// Copyright 2000-2022 CEA (www.cea.fr) IFPEN (www.ifpenergiesnouvelles.com)
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: Apache-2.0
//-----------------------------------------------------------------------------
/*---------------------------------------------------------------------------*/
/* Ptr.h                                                       (C) 2000-2006 */
/*                                                                           */
/* Classes diverses encapsulant des pointeurs.                               */
/*---------------------------------------------------------------------------*/
#ifndef ARCANE_UTILS_PTR_H
#define ARCANE_UTILS_PTR_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arcane/utils/ArcaneGlobal.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_BEGIN_NAMESPACE

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \ingroup Core
 * \brief Encapsulation d'un pointeur.
 *
 * Cette classe ne fait rien de particulier si ce n'est encapsulé un
 * pointeur d'un type quelconque. Elle sert de classe de base à d'autres
 * classes qui fournissent des fonctionnalités plus évoluées comme AutoRefT.
 *
 * Afin d'éviter des copies malencontreuses, le constructeur de copie et
 * les opérateurs de copie sont protégés.
 *
 * En mode débug, vérifie qu'on accède pas à un pointeur nul.
 *
 * Le paramètre template n'a pas besoin d'être défini. Cette classe peut donc
 * être instanciée pour un type opaque.
 *
 \since 0.4.39
 \date 16/07/2001
 */
template<class T>
class PtrT
{
 protected:

  //! Opérateur de copie
  const PtrT<T>& operator=(const PtrT<T>& from)
    { m_value = from.m_value; return (*this); }
  
  template <typename T2>
  const PtrT<T>& operator=(const PtrT<T2>& from)
      { m_value = from.get(); return (*this); }


  //! Affecte à l'instance la value \a new_value
  const PtrT<T>& operator=(T* new_value)
    { m_value = new_value; return (*this); }

  //! Construit une référence référant \a from
  PtrT(const PtrT<T>& from) : m_value(from.m_value) {}

  //! Construit une référence référant \a from
  template <typename T2>
  PtrT(const PtrT<T2>& from) : m_value(from.m_value) {}


 public:

  //! Construit une instance sans référence
  PtrT() : m_value(0) {}

  //! Construit une instance référant \a t
  explicit PtrT(T* t) : m_value(t) {}

  virtual ~PtrT() {}

 public:
  
  //! Retourne l'objet référé par l'instance
  inline T* operator->() const
    {
#ifdef ARCANE_CHECK
      if (!m_value)
        arcaneNullPointerError();
#endif
      return m_value;
    }

  //! Retourne l'objet référé par l'instance
  inline T& operator*() const
    {
#ifdef ARCANE_CHECK
      if (!m_value)
        arcaneNullPointerError();
#endif
      return *m_value;
    }

  /*! \brief Retourne l'objet référé par l'instance
  *
  * \warning En général, il faut être prudent lorsqu'on utilise cette
   * fonction et ne pas conservé le pointeur retourné.
   */
  inline T* get() const
    {
      return m_value;
    }

  inline bool isNull() const {
    return (m_value == NULL);
  }

 protected:
  
  T* m_value; //!< Pointeur sur l'objet référencé
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Compare les objets référencés par \a v1 et \a v2
 * La comparaison se fait pointeur par pointeur.
 * \retval true s'ils sont égaux
 * \retval false sinon
 */
template<class T1,class T2> inline bool
operator==(const PtrT<T1>& v1,const PtrT<T2>& v2)
{
  return v1.get() == v2.get();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Compare les objets référencés par \a v1 et \a v2
 * La comparaison se fait pointeur par pointeur.
 * \retval false s'ils sont égaux
 * \retval true sinon
 */
template<class T1,class T2> inline bool
operator!=(const PtrT<T1>& v1,const PtrT<T2>& v2)
{
  return v1.get() != v2.get();
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

ARCANE_END_NAMESPACE

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif  

