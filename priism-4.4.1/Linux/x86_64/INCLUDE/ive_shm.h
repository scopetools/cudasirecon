#ifndef UCSF_MSG_IVE_SHM_H
#define UCSF_MSG_IVE_SHM_H

/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare types and functions for shared memory operations.
 */


#include "ive_conf.h"


struct IVEArenaImpl;
typedef struct IVEArenaImpl* IVEArena;
struct IVEArenaAttrImpl;
typedef struct IVEArenaAttrImpl* IVEArenaAttr;

/*
 * These entities are shared between processes so use types with a size
 * that is the same for different operating systems or compilation models.
 */
#ifdef IVECONF_SHM64
typedef uint64_t IVEArenaSize;
typedef uint64_t IVEPtrRep;
#else
typedef uint32_t IVEArenaSize;
typedef uint32_t IVEPtrRep;
#endif

#ifdef IVECONF_SHMMUTEX_USE_PTHREAD
#include <pthread.h>
typedef pthread_mutex_t IVEShmMutex;
#elif defined(IVECONF_SHMMUTEX_USE_SYSVSEM)
typedef struct { int sem_id; } IVEShmMutex;
#elif defined(IVECONF_SHMMUTEX_USE_FCNTL)
typedef struct { IVEPtrRep placeholder; } IVEShmMutex;
#else
#error "Do not know how to define IVEShmMutex."
#endif

#ifdef IVECONF_SHMHOLD_USE_SYSVSEM
typedef struct { int sem_id; } IVEShmHold;
#elif defined(IVECONF_SHMHOLD_USE_FCNTL)
typedef struct { IVEPtrRep placeholder; } IVEShmHold;
#else
#error "Do not know how to define IVEShmHold."
#endif

#ifdef IVECONF_SHMRWLK_USE_PTHREAD
#include <pthread.h>
typedef pthread_rwlock_t IVEShmRWLock;
#elif defined(IVECONF_SHMRWLK_USE_SYSVSEM)
typedef struct { int sem_id; } IVEShmRWLock;
#elif defined(IVECONF_SHMRWLK_USE_FCNTL)
typedef struct { IVEPtrRep placeholder; } IVEShmRWLock;
#else
#error "Do not know how to define IVEShmRWLock."
#endif

#ifdef IVECONF_SHMSEM_USE_POSIXUNNAMED
#include <semaphore.h>
typedef sem_t IVEShmSemaphore;
#elif defined (IVECONF_SHMSEM_USE_SYSVSEM)
typedef struct { int sem_id; } IVEShmSemaphore;
#elif defined (IVECONF_SHMSEM_USE_POSIXNAMED)
typedef struct { char name[40]; } IVEShmSemaphore;
#elif defined (IVECONF_SHMSEM_USE_CYGWIN_HACK)
typedef struct {
    char mutex_name[50]; char event_name[50]; unsigned int count;
} IVEShmSemaphore;
#else
#error "Do not know how to define IVEShmSemaphore."
#endif


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

IVEArenaAttr IVECreateArenaAttr(void);
void IVEDestroyArenaAttr(IVEArenaAttr attr);
void IVESetArenaAttrSize(IVEArenaSize sz, IVEArenaAttr attr);
void IVESetArenaAttrChunk(IVEArenaSize chunk_sz, IVEArenaAttr attr);
IVEArena IVECreateArena(
    const char* filename, IVEArenaAttr attr, int fail_if_exist, int* p_status
);
void IVEDestroyArena(IVEArena arena);
IVEArena IVEAttachArena(const char* filename, int* p_status);
void IVEDetachArena(IVEArena arena);

IVEPtrRep IVEGetArenaRoot(IVEArena arena);
void IVESetArenaRoot(IVEPtrRep value, IVEArena arena);
int IVECompareSwapArenaRoot(
    IVEPtrRep old_value, IVEPtrRep new_value, IVEArena arena
);

void* IVEShmAlloc(IVEArenaSize sz, IVEArena arena);
void* IVEShmCalloc(IVEArenaSize n_elem, IVEArenaSize el_sz, IVEArena arena);
void* IVEShmRealloc(void* ptr, IVEArenaSize sz, IVEArena arena);
void IVEShmDealloc(void* ptr, IVEArena arena);
IVEPtrRep IVEPtrToRep(void* ptr, IVEArena arena);
void* IVERepToPtr(IVEPtrRep rep, IVEArena arena);

void IVEShmInitMutex(IVEShmMutex* p_mutex, int* p_status);
void IVEShmDestroyMutex(IVEShmMutex* p_mutex);
void IVEShmAcquireMutex(IVEShmMutex* p_mutex);
int IVEShmTryAcquireMutex(IVEShmMutex* p_mutex);
void IVEShmReleaseMutex(IVEShmMutex* p_mutex);

void IVEShmInitHold(IVEShmHold* p_hold, int* p_status);
void IVEShmDestroyHold(IVEShmHold* p_hold);
void IVEShmAcquireHold(IVEShmHold* p_hold);
int IVEShmTryAcquireHold(IVEShmHold* p_hold);
void IVEShmReleaseHold(IVEShmHold* p_hold);

void IVEShmInitRWLock(IVEShmRWLock* p_lock, int* p_status);
void IVEShmDestroyRWLock(IVEShmRWLock* p_lock);
void IVEShmAcquireReadLock(IVEShmRWLock* p_lock);
int IVEShmTryAcquireReadLock(IVEShmRWLock* p_lock);
void IVEShmReleaseReadLock(IVEShmRWLock* p_lock);
void IVEShmAcquireWriteLock(IVEShmRWLock* p_lock);
int IVEShmTryAcquireWriteLock(IVEShmRWLock* p_lock);
void IVEShmReleaseWriteLock(IVEShmRWLock* p_lock);

void IVEShmInitSemaphore(
    unsigned int initial_value, IVEShmSemaphore* p_sema, int* p_status
);
void IVEShmDestroySemaphore(IVEShmSemaphore* p_sema);
void IVEShmIncrementSemaphore(IVEShmSemaphore* p_sema);
void IVEShmDecrementSemaphore(IVEShmSemaphore* p_sema);
int IVEShmTryDecrementSemaphore(IVEShmSemaphore* p_sema);
int IVEShmGetSemaphoreValue(IVEShmSemaphore* p_sema);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
