#ifndef UCSF_MSG_FORTRAN_TYPES_H
#define UCSF_MSG_FORTRAN_TYPES_H

/* @(#) $Id: fortran_types.h,v 1.1.1.1 2001/07/12 23:12:59 eric Exp $ */
/* $Name:  $ */
/*
 * Define types for interfaces between Fortran and C/C++ code.
 */


#include "ive_conf.h"

#ifdef FTN_INT_TYPE
typedef FTN_INT_TYPE ftn_int_t;
#endif /* FTN_INT_TYPE */

#ifdef FTN_REAL_TYPE
typedef FTN_REAL_TYPE ftn_real_t;
typedef struct ftn_cmplx_t { ftn_real_t real, imag; } ftn_cmplx_t;
#endif /* FTN_REAL_TYPE */

#ifdef FTN_DOUBLE_TYPE
typedef FTN_DOUBLE_TYPE ftn_double_t;
typedef struct ftn_dcmplx_t { ftn_double_t real, imag; } ftn_dcmplx_t;
#endif /* FTN_DOUBLE_TYPE */

#ifdef FTN_CHAR_TYPE
typedef FTN_CHAR_TYPE ftn_char_t;
#endif /* FTN_CHAR_TYPE */

/*
 * The type used to represent the lengths of character arrays when passed
 * as arguments.
 */
#ifdef FTN_LENGTH_TYPE
typedef FTN_LENGTH_TYPE ftn_length_t;
#endif /* FTN_LENGTH_TYPE */

#ifdef FTN_LOGICAL_TYPE
typedef FTN_LOGICAL_TYPE ftn_logical_t;
#endif /* FTN_LOGICAL_TYPE */
#undef FTN_IS_FALSE
#ifdef FTN_FALSE_IS_ZERO
#define FTN_IS_FALSE(x) ((x) == 0)
#endif /* FTN_FALSE_IS_ZERO */
#ifdef FTN_FALSE_LOW_BIT_ZERO
#define FTN_IS_FALSE(x) (((x) & 1) == 0)
#endif /* FTN_FALSE_LOW_BIT_ZERO */

#endif /* include guard */
