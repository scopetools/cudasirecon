#ifndef UCSF_MSG_CONVERT_H
#define UCSF_MSG_CONVERT_H

#include "ive_conf.h"  /* uint32_t, uint16_t */

#define IMswab2(y) (((uint16_t)(y) >> 8) | (((uint16_t)(y) & 0xffu) << 8))
#define IMswab4(y) (((uint32_t)(y) >> 24) | (((uint32_t)(y) & 0xffu) << 24) | (((uint32_t)(y) & 0xff00u) << 8) | (((uint32_t)(y) & 0xff0000u) >> 8))

/*
 * Defines a prototype for converting to image data as seen by the user
 * (native byte ordering and possibly converted to floating-point) to the
 * representation stored in an image resource.  nc is the columns of data
 * to convert.  nr is the number of rows to convert.  nelem is the
 * number of elements per sample (two for complex data, one for others).
 */
typedef void (*IMToUserConverterInPlace)(
    void* inout, int nc, int nr, int nelem
);

/*
 * Similar to IMToUserConverterInPlace but for use where the operation can not
 * be done in place.  Takes two additional parameters.  ilskip is the number
 * of bits to skip at the beginning of each row from the input.  itskip is
 * the number of bits to skip at the end of each row from the input.  ilskip
 * and itskip are only non-zero for types that do not take a whole number of
 * bytes per element.  Assumes the input and output have a whole number of
 * bytes per row.
 */
typedef void (*IMToUserConverter)(
    void* out,
    const void* in,
    int nc,
    int nr,
    int nelem,
    int ilskip,
    int itskip
);

/*
 * Defines a prototype for converting from image data as seen by the user
 * to that stored in the image resource.  nc is the number of columns and nr
 * is the number of rows for the pixel array to be converted.  nelem is the
 * number of elements per pixel (two for complex data, one for others).  olskip
 * is the number of bits to skip at the beginning of each row from the output.
 * otskip is the number of bits to skip at the end of each row from the output.
 * olskip and otskip are only non-zero for types that do not take a whole
 * number of bytes per element.  Returns 0 if no overflow was detected;
 * otherwise returns one.
 */
typedef int (*IMFromUserConverter)(
    void* out,
    const void* in,
    int nc,
    int nr,
    int nelem,
    int olskip,
    int otskip
);

#ifdef __cplusplus
extern "C" {
#endif

void IMvswab4(uint32_t* inout, size_t count);
void IMvswab4_copy(uint32_t* out, const uint32_t* in, size_t count);
void IMvswab2(uint16_t* inout, size_t count);
void IMvswab2_copy(uint16_t* out, const uint16_t* in, size_t count);

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
