#ifndef UCSF_MSG_CONVERT_H
#define UCSF_MSG_CONVERT_H

/* @(#) $Id: convert.h,v 1.1 2002/01/18 00:24:52 eric Exp $ */
/* $Name:  $ */

#include "ive_conf.h"  /* uint32_t, uint16_t */

#define IMswab2(y) (((uint16_t)(y) >> 8) | ((uint16_t)(y) << 8))
#define IMswab4(y) (((uint32_t)(y) >> 24) | ((uint32_t)(y) << 24) | (((uint32_t)(y) & 0xff00u) << 8) | (((uint32_t)(y) & 0xff0000u) >> 8))

/*
 * Defines a prototype for converting to image data as seen by the user
 * (native byte ordering and possibly converted to floating-point) to the
 * representation stored in an image resource.  count is the number of pixels,
 * or if the data is complex, twice the number of pixels, to convert.
 */
typedef void (*IMToUserConverterInPlace)(void* inout, int count);

/*
 * Same as IMToUserConverter but for use where the operation can not be done
 * in place.
 */
typedef void (*IMToUserConverter)(void* out, const void* in, int count);

/*
 * Defines a prototype for converting from image data as seen by the user
 * to that stored in the image resource.  count is the number of pixels,
 * or if the data is complex, twice the number of pixels, to convert.
 * Returns 0 if no overflow was detected; otherwise one is returned.
 */
typedef int (*IMFromUserConverter)(void* out, const void* in, int count);

#ifdef __cplusplus
extern "C" {
#endif

void IMvswab4(uint32_t* inout, int count);
void IMvswab4_copy(uint32_t* out, const uint32_t* in, int count);
void IMvswab2(uint16_t* inout, int count);
void IMvswab2_copy(uint16_t* out, const uint16_t* in, int count);

/*
 * Returns a functions which is appropriate for converting data in the
 * resource representation (specified by mode, which describes the
 * representation of each pixel, and by swab which, if true, specifies,
 * that byte swapping should be done) to the user representation (which
 * is the native floating-point form if to_float is true).
 */
IMToUserConverterInPlace IMget_to_user_convert_ip(
    int mode, int to_float, int swab
);

IMToUserConverter IMget_to_user_convert(int mode, int to_float, int swab);

/*
 * Returns a function which is appropriate for converting data from the
 * user's representation to the resource representation.
 */
IMFromUserConverter IMget_from_user_convert(int mode, int to_float, int swab);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
