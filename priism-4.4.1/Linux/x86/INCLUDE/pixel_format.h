#ifndef UCSF_MSG_PIXEL_FORMAT_H
#define UCSF_MSG_PIXEL_FORMAT_H

typedef enum PixelFormat {
    PIXEL_INVALID,    /* possible return value for match */
    PIXEL_RGB_32,     /* 8 bits (unsigned) per channel, 32 bits wide */
    PIXEL_RGBA_32,    /* 8 bits (unsigned) per channel, 32 bits wide */
    PIXEL_UINT_8,
    PIXEL_INT_8,
    PIXEL_UINT_16,
    PIXEL_INT_16,
    PIXEL_UINT_32,
    PIXEL_INT_32,
    PIXEL_IEEE_32,
    PIXEL_IEEE_64,
    PIXEL_REAL_IMAG_INT_16,
    PIXEL_REAL_IMAG_IEEE_32,
    PIXEL_RGB_24     /* 8 bits (unisgned) per channel, 24 bit packed pixels */
} PixelFormat;

#endif /* include guard */
