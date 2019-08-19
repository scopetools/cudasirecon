#ifndef UCSF_MSG_IMAGEHEADER_H
#define UCSF_MSG_IMAGEHEADER_H

#include "ive_conf.h"
#include "IWApiConstants.h"
#include <stdio.h>


/* As unsigned, C0A0 (49312 decimal) */
#define DVID -16224
/* As unsigned, A0C0 (41152 decimal) */
#define DVID_SWAPPED -24384
#define IMAGE_EXTRA_COUNT 24
/*
 * At some point in the past the extra field shrunk.  It used to start
 * IMAGE_EXTRA_OFFSET bytes before its curent position.
 */
#define IMAGE_EXTRA_OFFSET 8
#define IMAGE_LABEL_COUNT 10
#define IMAGE_LABEL_SIZE 80
#ifdef BIGENDIAN
#define IMAGE_MACHINE_FORMAT 0
#elif defined(LITTLEENDIAN)
#define IMAGE_MACHINE_FORMAT 1
#else
#error "Do not know endianness"
#endif


typedef struct ImageHeader {
    int32_t ncrs[3];
    int32_t mode;
    int32_t ncrst[3];
    int32_t nxyz[3];
    float32_t cel[6];
    int32_t mapcrs[3];
    float32_t denmmm[3];
    int32_t nspg;
    int32_t next;
    union {
	struct {
	    int16_t ndvid;
	    int16_t nblank;
	} p;
	struct {
	    int8_t extra0[4];
	} e;
    } va;
    int32_t ntst;
    int8_t extra[IMAGE_EXTRA_COUNT];
    int16_t nint;
    int16_t nreal;
    int16_t nres;
    int16_t nzfact;
    union {
	struct {
	    float32_t min2;
	    float32_t max2;
	    float32_t min3;
	    float32_t max3;
	    float32_t min4;
	    float32_t max4;
	} p;
	struct {
	    int8_t extra2[24];
	} e;
    } vb;
    int16_t file_type;
    int16_t lens;
    int16_t n1;
    int16_t n2;
    int16_t v1;
    int16_t v2;
    union {
	struct {
	    float32_t min5;
	    float32_t max5;
	    int16_t num_times;
	    int16_t interleaved;
	} p;
	struct {
	    int8_t extra3[12];
	} e;
    } vc;
    float32_t tilt[3];
    union {
	struct {
	    int16_t num_waves;
	    int16_t wavelengths[IW_MAX_WAVE];
	    float32_t origz;
	    float32_t origxy[2];
	} p;
	struct {
	    float32_t xyz0[3];
	    int8_t map[4];
	    int8_t machst[4];
	    float32_t rms;
	} e;
    } vd;
    int32_t nlab;
    int8_t labels[IMAGE_LABEL_SIZE * IMAGE_LABEL_COUNT];
} ImageHeader;


#ifdef __cplusplus
extern "C" {
#endif


/*
 * Returns true if the header appears to be byte-swapped from the native
 * ordering.
 */
int IMcheck_header_swabbed(const ImageHeader* p_hdr);

/*
 * Return the size in bytes per pixel of a given mode.  Returns 0 if the
 * mode is invalid or uses a fractional number of bytes per sample.  This
 * function is deprecated and will be removed in a future release.
 */
int IMcompute_pixel_size(int mode);

/*
 * This is an extended version of IMcompute_pixel_size() which also handles
 * types which use a fractional number of pixels per sample.  A single sample
 * takes *p_num /( (float) *p_den) bytes.  *p_pad indicates how samples are
 * padded to a full number of bytes and will be IW_BIT_PAD_SAMPLE,
 * IW_BIT_PAD_ROW, or IW_BIT_PAD_SECTION.  Sets *p_num to zero, *p_den to one,
 * and *p_pad to IW_BIT_PAD_SAMPLE if mode is not one of the known pixel types.
 */
void IMcompute_pixel_size_ext(int mode, int* p_num, int* p_den, int* p_pad);

/*
 * Computes the x/y and z reductions factors for the resolution res and the
 * given z downsampling.  res is assumed to be nonnegative and z_downsampling
 * is assumed to be positive.
 */
void IMcompute_resolution_sizes(
    int res, int z_downsampling, int* xy_factor, int* z_factor
);

/*
 * Computes the pixel spacings for the x, y, and z directions.  A spacing
 * of zero is reported for a direction if the map value was invalid or
 * the sampling was zero.
 */
void IMcompute_spacing(
    const int32_t sampling[3],
    const int32_t map[3],
    const float32_t cell_length[3],
    float spacing[3]
);

/*
 * Resets all values in the header to zero excepct mapcrs[0-2] (set to
 * 1, 2, 3), cell[0-2] (set to 1.0f), cell[3-5] (set to 90.0f),
 * ndvid (set to DVID), nres (set to 1), nzfact (set to 1), num_times
 * (set to 1), num_waves (set to 1), and labels (filled with spaces).  Assumes
 * p_hdr is valid.
 */
void IMinitialize_header(ImageHeader* p_hdr);
/*
 * Exactly like IMinitialize_header, except ndvid, num_times, num_waves aren't
 * applicable and fills the map with 'MAP ' and machst with the little- or
 * big-endian stamp as appropriate.
 */
void IMinitialize_emheader(ImageHeader* p_hdr, int is_little_endian);

/*
 * Writes a human-readable version of most of the fields in the main header.
 * assumes p_hdr is valid.
 */
void IMprint_header(const ImageHeader* p_hdr, FILE* fo);
void IMprint_emheader(const ImageHeader* p_hdr, FILE* fo);

/*
 * Performs the appropriate byte swapping on all entries in the main header.
 */
void IMswab_header(ImageHeader* p_hdr);
void IMswab_emheader(ImageHeader* p_hdr);

/*
 * Returns IW_EXTHDR_PRIISM0 if the extended header should be interpreted as
 * a fixed number of integers and floating-point values per section.  Returns
 * IW_EXTHDR_CRYSTAL if the extended header should be interpreted as
 * crystallographic symmetry information.
 */
int IMget_exthdr_format(const ImageHeader* p_hdr, int header_format);

/*
 * Return zero if none of the dimensions (nx, ny, nz, nt, nw) changed during
 * the conversion; otherwise, return one.
 */
int IMconvert_hdr_em0_to_priism0(ImageHeader* p_hdr);
int IMconvert_hdr_priism0_to_em0(ImageHeader* p_hdr, int stored_byte_order);

void IMget_zwt(
    const ImageHeader* p_hdr, 
    int nsec,
    int header_format,
    int* p_nz,
    int* p_nw,
    int* p_nt,
    int* p_intlv
);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
