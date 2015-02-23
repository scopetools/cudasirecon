#ifndef UCSF_MSG_IVE_THREAD_H
#define UCSF_MSG_IVE_THREAD_H

/* @(#) $Id: ive_thread.h,v 1.5 2008/04/07 22:29:18 eric Exp $ */
/* $Name:  $ */
/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare types and functions for making the library thread-safe.
 * These functions do nothing unless you link against the thread-safe
 * library and define IVE_MTSAFE.
 */


#include "ive_conf.h"


#ifdef IVE_MTSAFE
#ifdef IVECONF_PTHREAD
#include <pthread.h>
#define IVE_ONCE_INIT PTHREAD_ONCE_INIT
typedef pthread_once_t IVEThreadOnceFlag;
typedef pthread_key_t IVEThreadKey;
typedef pthread_mutex_t IVEThreadMutex;
#elif defined(IVECONF_WIN32)
#define IVE_ONCE_INIT (0)
typedef PVOID IVEThreadOnceFlag;
typedef DWORD IVEThreadKey;
typedef CRITICAL_SECTION IVEThreadMutex;
#else
#error "Do not know how to define IVEThread types"
#endif /* else clause IVECONF_PTHREAD */

#else
#define IVE_ONCE_INIT {0}
typedef struct { int set; } IVEThreadOnceFlag;
struct IVEThreadKeyImpl;
typedef struct IVEThreadKeyImpl* IVEThreadKey;
typedef struct { int dummy; } IVEThreadMutex;

#endif /* else clause IVE_MTSAFE */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int IVEThreadCountProcessors(void);
void IVEThreadCallOnce(void (*init_routine)(void), IVEThreadOnceFlag* p_once);
void IVEThreadInitKey(IVEThreadKey* p_key, int* p_status);
void IVEThreadDestroyKey(IVEThreadKey key);
void IVEThreadSetKeyValue(void* value, IVEThreadKey key, int* p_status);
void* IVEThreadGetKeyValue(IVEThreadKey key);
void IVEThreadInitMutex(IVEThreadMutex* p_mutex, int* p_status);
void IVEThreadDestroyMutex(IVEThreadMutex* p_mutex);
void IVEThreadAcquireMutex(IVEThreadMutex* p_mutex);
int IVEThreadTryAcquireMutex(IVEThreadMutex* p_mutex);
void IVEThreadReleaseMutex(IVEThreadMutex* p_mutex);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
