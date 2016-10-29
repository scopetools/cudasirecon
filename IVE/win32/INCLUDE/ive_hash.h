#ifndef UCSF_MSG_IVE_HASH_H
#define UCSF_MSG_IVE_HASH_H

/* @(#) $Id: ive_hash.h,v 1.2 2007/09/24 17:48:38 eric Exp $ */
/* $Name:  $ */
/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */
/*
 * Declare a utility function for hashing arbitrary data to an integer
 * value.  hash_size (the possible range of values for the hash is 0 to
 * hash_size - 1) is restricted to be a power of 2.
 */


#include "ive_conf.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int IVEHash(const uint8_t* data, uint32_t length, int hash_size);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
