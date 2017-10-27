#ifndef UCSF_MSG_IVE_PROCESSBASE_H
#define UCSF_MSG_IVE_PROCESSBASE_H

/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare functions for use in writing Fortran wrappers for C.
 */

#include <stddef.h>  /* size_t */

typedef void (*IVEFatalHandler)(const char* message, size_t length);
typedef void (*IVEWarningHandler)(const char* message, size_t length);


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void IVERegisterCleanupRoutine(
    void (*func)(void*), void* param, int* p_status
);

void IVEReportFatal(const char* message);

IVEFatalHandler IVESetFatalHandler(IVEFatalHandler handler);

IVEFatalHandler IVEGetFatalHandler(void);

void IVEReportWarning(const char* message);

IVEWarningHandler IVESetWarningHandler(IVEWarningHandler handler);

IVEWarningHandler IVEGetWarningHandler(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
