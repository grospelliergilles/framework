﻿// -*- tab-width: 2; indent-tabs-mode: nil; coding: utf-8-with-signature -*-
#include <gtest/gtest.h>

#include "arccore/base/ReferenceCounter.h"
#include "arccore/base/Ref.h"

#include <iostream>

using namespace Arccore;

namespace
{
struct StatInfo
{
  bool is_destroyed = false;
  int nb_add = 0;
  int nb_remove = 0;
  // Vérifie que 'is_destroyed' est vrai et qu'on a fait
  // le bon nombre d'appel à 'nb_add' et 'nb_remove'.
  bool checkValid(int nb_call)
  {
    if (nb_add!=nb_call){
      std::cout << "Bad nb_add";
      return false;
    }
    if (nb_remove!=nb_call){
      std::cout << "Bad nb_remove";
      return false;
    }
    return is_destroyed;
  }
};
}

class Simple1
{
 public:
  Simple1(StatInfo* stat_info) : m_nb_ref(0), m_stat_info(stat_info){}
  ~Simple1(){ m_stat_info->is_destroyed = true; }
 public:
  void addReference()
  {
    ++m_nb_ref;
    ++m_stat_info->nb_add;
    std::cout << "ADD REFERENCE r=" << m_nb_ref << "\n";
  }
  void removeReference()
  {
    --m_nb_ref;
    ++m_stat_info->nb_remove;
    std::cout << "REMOVE REFERENCE r=" << m_nb_ref << "\n";
    if (m_nb_ref==0){
      std::cout << "DESTROY!\n";
      delete this;
    }
  }
  Int32 m_nb_ref;
  StatInfo* m_stat_info;
};

namespace Arccore
{
template<>
class RefTraits<Simple1>
{
 public:
  typedef ReferenceCounterWrapper<Simple1> ImplType;
};
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

template<typename RefType> void
_doTest1(const RefType& ref_type)
{
  {
    RefType s3;
    {
      RefType s1(ref_type);
      RefType s2 = s1;
      {
        s3 = s2;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

// Teste si le compteur de référence détruit bien l'instance.
TEST(ReferenceCounter, Misc)
{
  typedef ReferenceCounter<Simple1> Simple1Reference;
  {
    StatInfo stat_info;
    _doTest1(Simple1Reference(new Simple1(&stat_info)));
    ASSERT_TRUE(stat_info.checkValid(4)) << "Bad destroy1";
  }
}

// Teste si le compteur de référence détruit bien l'instance.
TEST(ReferenceCounter, Ref)
{
  typedef Ref<Simple1> RefSimple1;
  {
    StatInfo stat_info;
    _doTest1(makeRef(new Simple1(&stat_info)));
    ASSERT_TRUE(stat_info.checkValid(4)) << "Bad destroy2";
  }
}
