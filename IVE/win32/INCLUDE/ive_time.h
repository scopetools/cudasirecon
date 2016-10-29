#ifndef UCSF_MSG_IVE_TIME_H
#define UCSF_MSG_IVE_TIME_H

/* @(#) $Id: ive_time.h,v 1.5 2012/06/02 03:04:36 eric Exp $ */
/* $Name:  $ */
/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare types and functions for timing.
 */

#include "ive_conf.h"


#undef IVE_TIME_INTERVAL_MIN
#undef IVE_TIME_INTERVAL_MAX


/*
 * IVETimeValue is the type used to store a representation of the current time.
 * IVETimerInterval is the signed integral type used to store time intervals;
 * it must have at least 32 bits.
 */
#if defined(IVECONF_XOPEN) || defined(IVECONF_BSD)
#include <sys/time.h>
#include <limits.h>
typedef struct timeval IVETimeValue;
typedef int IVETimeInterval;
#define IVE_TIME_INTERVAL_MIN INT_MIN
#define IVE_TIME_INTERVAL_MAX INT_MAX
#elif defined(IVECONF_WIN32)
#include <windows.h>
#include <limits.h>
typedef FILETIME IVETimeValue;
typedef int IVETimeInterval;
#define IVE_TIME_INTERVAL_MIN INT_MIN
#define IVE_TIME_INTERVAL_MAX INT_MAX
#else
#error "Do not know how to define IVETimeValue and IVETimerInterval"
#endif


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void IVESleep(IVETimeInterval intvl_us);
void IVEGetCurrTime(IVETimeValue* p_time);
IVETimeInterval IVECalcTimeDiff(
    const IVETimeValue* p_t2, const IVETimeValue* p_t1, int* p_status
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
