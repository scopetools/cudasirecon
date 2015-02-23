#ifndef UCSF_MSG_IVE_FORTRAN_H
#define UCSF_MSG_IVE_FORTRAN_H

/* @(#) $Id: ive_fortran.h,v 1.3 2007/09/24 17:48:38 eric Exp $ */
/* $Name:  $ */
/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare functions for use in writing Fortran wrappers for C.
 */


#include "fortran_types.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 * Result should be freed with free.
 */
char* IVEFtnStringToCString(const ftn_char_t* fstr, ftn_length_t flength);

void IVECStringToFtnString(
    const char* cstr, ftn_length_t flength, ftn_char_t* fstr
);

ftn_int_t IVEPtrToFtnId(void* ptr);

void* IVEFtnIdToPtr(ftn_int_t id);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
