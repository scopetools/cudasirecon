#ifndef UCSF_MSG_IVE_ERROR_H
#define UCSF_MSG_IVE_ERROR_H

/* @(#) $Id: ive_error.h,v 1.3 2007/09/24 17:48:36 eric Exp $ */
/* $Name:  $ */
/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declares error constants and error handling functions.
 */


enum {
    IVE_ERR_NONE = 0

    /* Runtime conditions */
    , IVE_ERR_NO_MEM = 1
    , IVE_ERR_NO_SPC = 2
    , IVE_ERR_NO_RESOURCES = 3
    , IVE_ERR_IN_USE = 10
    , IVE_ERR_TIMEOUT = 11
    , IVE_ERR_IO = 20
    , IVE_ERR_NO_PERMISSION = 30
    , IVE_ERR_EXIST = 31
    , IVE_ERR_OVERFLOW = 40
    , IVE_ERR_UNSUPPORTED = 50

    /* Problems with passed arguments */
    , IVE_ERR_INVALID = 512
    , IVE_ERR_BAD_FILENAME = 513

    /* Miscellaneous */
    , IVE_ERR_CORRUPT = 1024
    , IVE_ERR_INTERNAL = 1025
    , IVE_ERR_UNRECOGNIZED = 1280
};

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void IVEErrorToString(int code, int max_length, char* str, int* p_full_length);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
