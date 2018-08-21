/*---------------------------------------------------------------------------*/
/* IThreadImplementation.h                                     (C) 2000-2009 */
/*                                                                           */
/* Interface d'un service implémentant le support des threads.               */
/*---------------------------------------------------------------------------*/
#ifndef ARCCORE_CONCURRENCY_ITHREADIMPLEMENTATION_H
#define ARCCORE_CONCURRENCY_ITHREADIMPLEMENTATION_H
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#include "arccore/concurrency/ConcurrencyGlobal.h"
#include "arccore/base/BaseTypes.h"

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

namespace Arccore
{

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*!
 * \brief Interface d'un service implémentant le support des threads.
 */
class ARCCORE_CONCURRENCY_EXPORT IThreadImplementation
{
 public:
  virtual ~IThreadImplementation(){}
 public:
  virtual void initialize() =0;
 public:
  virtual ThreadImpl* createThread(IFunctor* f) =0;
  virtual void joinThread(ThreadImpl* t) =0;
  virtual void destroyThread(ThreadImpl* t) =0;

  virtual void createSpinLock(Int64* spin_lock_addr) =0;
  virtual void lockSpinLock(Int64* spin_lock_addr,Int64* scoped_spin_lock_addr) =0;
  virtual void unlockSpinLock(Int64* spin_lock_addr,Int64* scoped_spin_lock_addr) =0;

  virtual MutexImpl* createMutex() =0;
  virtual void destroyMutex(MutexImpl*) =0;
  virtual void lockMutex(MutexImpl* mutex) =0;
  virtual void unlockMutex(MutexImpl* mutex) =0;

  virtual Int64 currentThread() =0;

  virtual IThreadBarrier* createBarrier() =0;
};

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

} // End namespace Arccore

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

#endif
