#ifndef UCSF_MSG_IWAPICONSTANTS_H
#define UCSF_MSG_IWAPICONSTANTS_H

/* @(#) $Id: IWApiConstants.h,v 1.9 2009/11/10 01:06:52 eric Exp $ */
/* $Name:  $ */
/*-----------------------------------------------------------------------------
    Copyright (C) 1993
    Macromolecular Structure Group of Biochemistry Dept. at University of
    California at San Francisco.
    These coded instructions, statements, and computer programs comprise
    unpublished proprietary information of the Macromolecular Structure Group of
    the Biochemistry Department at University of California at San Francisco, 
    and are protected by Federal copyright law.  They may not be disclosed
    to third parties, copied or duplicated in any form - in whole or in part -
    without the prior written consent of Macromolecular Structure Group of
    Biochemistry Department at University of California at San Francisco.
-------------------------------------------------------------------------*/

#define  IW_SUCCESS             1
#define  IW_FAILURE             0
#define  IW_ERROR               -1

/* Constants for designating bank channels. */
#define IW_EMPTY                 -1
#define IW_ALL_WINDOWS           -1
#define IW_ALL_BANKS             -1
#define IW_MAX_BANK              5
#define IW_ALL_WAVES             -1
#define IW_WAVE1                 0
#define IW_WAVE2                 1
#define IW_WAVE3                 2
#define IW_WAVE4                 3
#define IW_WAVE5                 4
#define IW_MAX_WAVE              5

/* for IMWLapply_to_wave */
#define IW_APPLY_SYNC            0
#define IW_APPLY_NONE            1
#define IW_APPLY_ALL             2

/* Constants for specifing color */
#define IW_PSEUDO                1
#define IW_TRUE_COLOR            0
#define IW_COLOR_NAME_LEN        12
#define IW_NUM_TRUE_COLOR        7

#define IW_BLACK                 0
#define IW_RED                   1
#define IW_GREEN                 2
#define IW_YELLOW                (IW_RED + IW_GREEN)
#define IW_BLUE                  4
#define IW_PINK                  (IW_RED + IW_BLUE)
#define IW_CYAN                  (IW_GREEN + IW_BLUE)
#define IW_WHITE                 (IW_RED + IW_GREEN + IW_BLUE)
#define IW_SANDY_BROWN           8
#define IW_BROWN                 9
#define IW_LT_GREY               10
#define IW_PURPLE                11
#define IW_INVISIBLE             12
#define IW_NUM_COLORS            12

/* Data Types. */
#define  IW_AS_IS                -1
#define  IW_BYTE                  0
#define  IW_SHORT                 1
#define  IW_FLOAT                 2
#define  IW_COMPLEX_SHORT         3
#define  IW_COMPLEX               4
#define  IW_EMTOM                 5
#define  IW_USHORT                6
#define  IW_LONG                  7

/*
 * Possible values for the file type (see IMAlDat and IMRtDat).  These
 * are new in IVE 4; previous versions referred to them directly by number.
 * IM_MULTIPOSITION was added in IVE 4.1.0; IM_PUPIL_FUNCTION was added in
 * IVE 4.1.8.
 */
#define IM_NORMAL_IMAGES          0
#define IM_TILT_SERIES            1
#define IM_STEREO_TILT_SERIES     2
#define IM_AVERAGED_IMAGES        3
#define IM_AVERAGED_STEREO_PAIRS  4
#define IM_EM_TILT_SERIES         5
#define IM_MULTIPOSITION          20
#define IM_PUPIL_FUNCTION         8000

/* Monitor Matrix Display Layout Format */
#define TEXT_FMT        1
#define DATA_FMT        0

/*
   Monitor window decoration style (for IWAlDisWinBdr, IWRtDisWinBdr).
   Full border includes the title bar and the resize handles.  The resize
   border only has the resize handles.  Windows without borders have
   no window manager decorations.
*/
#define IW_MON_FULL_BRDR 0
#define IW_MON_RESIZE_BRDR 1
#define IW_MON_NO_BRDR   2

/*
   Monitor menu/toolbar display style (for IWAlDisWinTool, IWRtDisWinTool).
   No menus or tools are shown with IW_MON_NO_MENU_TOOLS.  The other constants
   can be combined as a bitwise or to flag what will be shown.
*/
#define IW_MON_NO_MENU_TOOLS   0
#define IW_MON_SHOW_MENU       1
#define IW_MON_SHOW_LEFTTOOLS  2
#define IW_MON_ALL_MENU_TOOLS  3

#define  IW_MIN_WIN_SIZE          128
#define  IW_MAX_WIN_SIZE          1024


/*
 * Specify which buffers to read out for IWCaptureImage.
 */
#define IW_LEFT_IMAGE    1
#define IW_RIGHT_IMAGE   2
#define IW_NORMAL_IMAGE  3

/*
 * Specify which format to use for each pixel in the result of IWCaptureImage.
 */
#define IW_RGBA_PIXEL      0
#define IW_LUMINANCE_PIXEL 1


typedef struct IW_MRC_Header {
  int nx,ny,nz;               /* nz : nfocal*nwave*ntime */
  int mode;
  int nxst, nyst, nzst;
  int mx, my, mz;
  float xlen, ylen, zlen;
  float alpha, beta, gamma;
  int mapc, mapr, maps;
  float amin, amax, amean;
  int ispg, inbsym;

  short nDVID,nblank;      /* CB: nDVID value. nblank preserves byte boundary */
  int  ntst;               /* Hans: add ntst to record time domain offset */
  char ibyte[24];
  short nint,nreal;
  short nres,nzfact;
  float min2,max2,min3,max3,min4,max4;
  short file_type, lens, n1, n2, v1, v2; /* From daa2:[agard.image]imsubs.for*/
      /*file_type = idtype*/
  float min5,max5;
  short num_times;
  short interleaved;   /* Zero => not interleaved. That means z changes
                            fastest, then time, then waves;
                          1 => interleaved. That means wave changes fastest,
                            then z, then time. */
  float tilt_x, tilt_y, tilt_z;
  short num_waves, iwav1, iwav2, iwav3, iwav4, iwav5;
  float zorig, xorig, yorig;
  int nlab;
  char label[800];
} IW_MRC_HEADER, *IW_MRC_HEADER_PTR;

typedef struct IW_complex {
    float real;
    float imaginary;
} IW_COMPLEX_VALUE, *IW_COMPLEX_PTR;

/* for IMWLcomplex_disp */
#define IW_COMPLEX_AMPLITUDE     0x00010000
#define IW_COMPLEX_REAL          0x00020000
#define IW_COMPLEX_IMAGINARY     0x00030000
#define IW_COMPLEX_PHASE         0x00040000

/*
 * Used to define accuracy for scaling done in IWScaleImage.
 * IW_SCALING_DONTCARE is the default - exact for byte and short
 * types but increasingly inaccurate for floating point data as
 * the exponent gets farther away from 1.  IW_SCALING_EXACT is
 * exact for all types (but generally involves a performance penalty).
 */
#define IW_SCALING_DONTCARE 0x00000000
#define IW_SCALING_EXACT    0x01000000

typedef unsigned char      *IW_BYTE_MEM_PTR;
typedef short              *IW_SHORT_MEM_PTR;
typedef float              *IW_FLOAT_MEM_PTR;

#define IW_SAVE_BUF                 (1 << 0)
#define IW_BUF_CREATED              (1 << 1)
#define IW_IMG_CREATED              (1 << 2)

#define IW_WIN_CREATED               680     /* MAGIC NUMBER */
#define IW_WIN_KILLED                681     /* MAGIC NUMBER */
#define IW_BANK_DISPLAYED            689     /* MAGIC NUMBER */
#define IW_BANK_DISPLAY_CHANGED      690     /* MAGIC NUMBER */
#define IW_BANK_KILLED               691     /* MAGIC NUMBER */
#define IW_WINDOW_GEO_CHANGED        692     /* MAGIC NUMBER */
#define IW_WIN_GR_KILLED             693     /* MAGIC NUMBER */
#define IW_WAVE_GR_KILLED            694     /* MAGIC NUMBER */

/* Macro Function Names */
#define IW_MAX(x,y)          ((x > y) ? x : y) 
#define IW_MIN(x,y)          ((x > y) ? y : x) 
#define IW_ABS(x)            ((x > 0) ? x : -x)

typedef float IW_VECTOR[3];
typedef float *IW_VECTOR_PTR;
typedef IW_VECTOR IW_MATRIX[3];
typedef IW_VECTOR *IW_MATRIX_PTR;

typedef int IW_INT_VECTOR[3];
typedef int *IW_INT_VECTOR_PTR;
typedef IW_INT_VECTOR IW_INT_MATRIX[3];
typedef IW_INT_VECTOR *IW_INT_MATRIX_PTR;

/*
 * Define indices into the ClientMessage data.s array for Monitor display
 * change events.  Use 2 and 5 for backward compatibility.
 */
#define IW_CM_WIN_ID              2
#define IW_CM_CHANGE_OR_SYN       5 /* type of change */

#define GW_ATTEND                 0
#define GW_ATTACH                 1
#define GW_DEATTACH               2
#define GW_DELETEGR               3

/* Extended header definitions -- new additions */

#define UNKNOWN_EXT_HEADER_TYPE   0
#define UCSF_EXT_HEADER_TYPE1     1
#define API_EXT_HEADER_TYPE1      2
#define API_EXT_HEADER_TYPE2      3
#define UCSF_EXT_HEADER_TYPE2     4

#define UCSF_EXT_HDR_1_NINTS      1
#define UCSF_EXT_HDR_1_NFLOATS    1
#define UCSF_EXT_HDR_2_NINTS      3
#define UCSF_EXT_HDR_2_NFLOATS    1
#define API_EXT_HDR_1_NINTS       0
#define API_EXT_HDR_1_NFLOATS     2
#define API_EXT_HDR_2_NINTS       8
#define API_EXT_HDR_2_NFLOATS     32

/* NB: use of these values as array indices for the floating-point  */
/*     part of the extended header is a temporary KLUGE. */
/*     WHEN APPLYING THIS KLUGE IN FORTRAN, 1 MUST BE ADDED TO VALUE IN CODE! */
#define PHOTOSENSOR_READING 0
#define TIME_STAMP_SECS     1
#define STAGE_X_COORD       2
#define STAGE_Y_COORD       3
#define STAGE_Z_COORD       4
#define MIN_INTEN           5
#define MAX_INTEN           6
#define MEAN_INTEN          7
#define EXP_TIME            8
#define ND_FILTER           9
#define EX_WAVELEN          10
#define EM_WAVELEN          11
#define INTEN_SCALING       12
#define ENERGY_CONV_FACTOR  13
#define PHOTOSENSOR_NCONV   14

/* array indices for extended header integer values */
#define PHOTOMETRIC_INTERP_INT    0  /*0=white is zero,1=black is zero,2=RGB,3=palette,4=transparency (same as TIFF)*/

/* code values for wavelength */
#define RATIO_IMAGE_WAVELENGTH         10
#define PRODUCT_IMAGE_WAVELENGTH       20
#define ADDITION_IMAGE_WAVELENGTH      30
#define SUBTRACTION_IMAGE_WAVELENGTH   40
#define DIC_IMAGE_WAVELENGTH           50
#define PHASE_IMAGE_WAVELENGTH         60
#define BLEND_IMAGE_WAVELENGTH         70


/* image sequence definitions */
#define ZTW_SEQUENCE   0    /* "non-interleaved", by defintion       */
#define WZT_SEQUENCE   1    /* "interleaved", from R3D and others    */
#define ZWT_SEQUENCE   2    /* new sequence. Unsupported as of 11/97 */

#endif /* include guard */
