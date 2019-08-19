#ifndef UCSF_MSG_IWCONSTANTS_H
#define UCSF_MSG_IWCONSTANTS_H

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

/*
 * Files including this file should include ive_standards.h before
 * any headers (to appropriately set feature-set macros for IVEShmMutex
 * and IVEShmSemaphore).
 */
#include "ive_shm.h"    /* IVEPtrRep, IVEShmMutex, IVEShmSemaphore */
#include <X11/Xlib.h>   /* Colormap, Window, XID */
#include <sys/types.h>  /* pid_t */

/*
 * Allow compatibility checking for statically linked executables.
 * See compat_major and compat_minor fields in IW_shared_mem_struct and
 * IWAttach and IWStart in IWStart.c.
 */
#define IW_COMPAT_MAJOR           4
#define IW_COMPAT_MINOR           4

#define IW_DEFAULT_SHM_SIZE       300
#define IW_DEFAULT_WORKING_UNIT   16
#define IW_DEFAULT_WORKING_SET    96
#define IW_DEFAULT_WIN_SIZE       8
#define IW_DEFAULT_SEC_SIZE       1000
#define IW_MIN_SHM_SIZE           30
#ifdef IVECONF_SHM64
/* IRIX64 has a user space limit of 2^40, others are probably more generous. */
#define IW_MAX_SHM_SIZE           1048575
#else
#define IW_MAX_SHM_SIZE           4095
#endif
#define IW_DEFAULT_PGTHRESH       0

#define IW_MAX_COLOR_VALUE         (1<<16)
#define IW_NUM_PSEUDO_LUT          0
#define IW_NUN_TRUE_LUT            0

/****************************************************************************
  GLOBAL CONSTANTS FOR THE IAWS LIBRARY ROUTINES. 
 ****************************************************************************/

/* Some memonic constants. */
#define  IW_TRUE                1
#define  IW_FALSE               0
#define  IW_KEEP                0xFFFFFFFF

#define  IW_NO_ERROR            0 
#define  IW_UP                  1
#define  IW_DOWN                0

#define IW_FILE_NAME_SIZE       100
#define IW_TITLE_SIZE           100

#define IW_WIN_WIDTH            512
#define IW_WIN_HEIGHT           512 
#define IW_XYZ_DIM              3
#define IW_ZTW                  0  /* ZTW. Set wave, take all focus points,
                                      next wave. Waves non-interleaved. */
#define IW_WZT                  1  /* WZT. Set focus, take all waves,
                                      change focus. Waves interleaved. */
#define IW_ZWT                  2  /* */

#define IW_BANK1                 0
#define IW_BANK2                 1
#define IW_BANK3                 2
#define IW_BANK4                 3
#define IW_BANK5                 4

#define PrevImage                0
#define NextImage                1
#define FirstImage               2
#define LastImage                3


#define IW_MAX_ZOOM              32

#define IW_TOP                   0
#define IW_BOTTOM                -1
#define IW_CURRENT               -2
#define IW_FIND_WID              -3
#define IW_DISPLAY_BUFFER        -4

#define  IW_BY_SECTION            0
#define  IW_BY_Z                  1
#define  IW_BY_TIME               2
#define  IW_BY_WAVE               3

/* Hardware configurations. */
#define IW_SGL_RGB_BUFFER         0
#define IW_DBL_RGB_BUFFER         1
#define IW_Z_RGB_BUFFER           2
#define IW_NO_RGB_BUFFER          3
#define IW_NO_RGB_OVERLAY         4

/*
 * Define indices into ClientMessage data.l array for communication with
 * image windows.
 */
#define IW_CM_ACTION              0
#define IW_CM_ARG1                1
#define IW_CM_ARG2                2
#define IW_CM_ARG3                3
#define IW_CM_ARG4                4

/* Monitor Constants. */
#define   IW_MON_SYNC             0
#define   IW_MON_NO_SYNC          1

#define SWAP                      1
#define NO_SWAP                   0

/*
 * Define allowed values for data.l[IW_CM_ACTION] in the ClientMessage event.
 */
#define  IW_MON_PAINT_WIN         1    /* no arguments */
#define  IW_MON_SET_OVERLAY       2    /* one argument:  color */
#define  IW_MON_CREATE_WIN        3    /* one argument:  window number */
#define  IW_MON_CAPTURE_IMAGE     4    /* two arguments:  offset low and
                                        * high */

#define  IW_MON_CREATE_SEC        5    /* no arguments */
#define  IW_MON_CLEAR_WIN         6    /* one argument:  0 clear image;
                                        * 1 clear overlay */
#define  IW_MON_UPDATE_TITLE      7    /* no arguments */
#define  IW_MON_AL_SIZE           8    /* one argument:  number of sections */
#define  IW_MON_DELETE_WIN        9    /* one argument:  0 overwrite with
                                        * non-scratch: 1 or 3 delete
					* completely; 2 overwrite with
					* scratch  */

#define  IW_MON_BORDER_STYLE     10    /* one argument: 0 change to window
					* manager borders; 1 change to
					* toolbars */
#define  IW_MON_MOVE_RESIZE      11    /* no arguments */
#define  IW_MON_RAISE            12    /* no arguments */
#define  IW_MON_LOWER            13    /* no arguments */
#define  IW_MON_PREV_SEC         14    /* no arguments */

#define  IW_MON_NEXT_SEC         15    /* no arguments */
#define  IW_MON_LAST_SEC         16    /* no arguments */
#define  IW_MON_FIRST_SEC        17    /* no arguments */
#define  IW_MON_CHANGE_DIS_MODE  18    /* no arguments */
#define  IW_MON_ICONIFY_CHANGE   19    /* no arguments */

#define  IW_MON_DRAW_FAST_WINDOW_GRAPHICS       20    /* no arguments */
#define  IW_MON_PREPARE_FAST_WINDOW_GRAPHICS    21    /* no arguments */
#define  IW_MON_END_FAST_WINDOW_GRAPHICS        22    /* no arguments */
#define  IW_UPDATE_5D_MODEL      23    /* two arguments:  offset low and
                                        * high */
#define  IW_MON_CHANGE_OV_COLORMAP              24    /* no arguments */

#define  IW_MON_CHANGE_COLORMAP  25    /* no arguments */
#define  IW_MON_REG_EVENT        26    /* three arguments: destination window,
                                        * event mask, buttons */
#define  IW_MON_UN_REG_EVENT     27    /* three arguments: destination window,
                                        * event mask, buttons */
#define  IW_MON_REG_DISPLAY_CHANGE              28    /* one argument:
                                                       * destination window */
#define  IW_MON_UN_REG_DISPLAY_CHANGE           29    /* one argument:
                                                       * destination window */

#define  IW_MON_DIS_GR_GID       30    /* one argument:  graphic id */
#define  IW_MON_DIS_LIST         31    /* three arguments:  number of graphics,
                                        * offset low and high */
#define  IW_MON_SET_BUF_IMG      32    /* one argument:  on/off */
#define  IW_MON_CHANGE_STEP      33    /* no arguments */
#define  IW_MON_MOVE_POINTER     34    /* two arguments:  x, y */
#define  IW_MON_CHANGE_CURSOR    35    /* one argument: shape */
#define  IW_MON_CHANGE_WAVELENGTH               36    /* no arguments */
#define  IW_MON_CHANGE_RESOLUTION               37    /* one argument:
							 resolution level */
#define  IW_MON_NUM_FUNC                        38

#define  IW_ANGSTROMS             0
#define  IW_MICRONS               1

#define  IW_BOGUS_MIN             10000000.0
#define  IW_BOGUS_MAX             -10000000.0

typedef struct memStruct_ {
    IVEPtrRep prev_off;         /* offset to previous memStruct in list */
    IVEPtrRep next_off;         /* offset to next memStruct in list */
    IVEPtrRep mem_off;          /* offset of data array */
    long mem_size;              /* size of data array in bytes */
    int mem_lock;               /* flag for locking data in SHM */
    int type;                   /* data or display image */
    int wid;                    /* associated window ID */
    int sec;                    /* sec # of the associated window */
    int res;                    /* resolution level */
} memStruct;

typedef struct memQueue_ {
    IVEPtrRep Head;         /* offset of first member in the queue */
    IVEPtrRep Last;       /* offset of last member in the queue */
} memQueue;                 /* queue is a double link list */

#define IW_MEM_NO_FREE_IN_NO_MATCH   0
#define IW_MEM_FREE_IN_NO_MATCH      1
/*************************************************************************
  STRUCTURE DEFINITIONS FOR LINKED LISTS
 *************************************************************************/

typedef struct IWL_struct  *IWL_STRUCT_PTR; /* Overall structure for shared 
                                               memory. */
typedef struct IW_WL_  *IW_WL_PTR;  /* Window List. */
typedef struct IW_BAL_ *IW_BAL_PTR; /* Bank Attribute List. */
typedef struct IW_IAL_ *IW_IAL_PTR; /* Image Attribute List. This is the same
                                       as IW_IMAGE_INFO_PTR */
typedef struct IW_GIL_ *IW_GIL_PTR; /* Geometry Info List. This is the same as
                                       IW_IMAGES_INFO_PTR */

typedef char               IW_FILE_STR[IW_FILE_NAME_SIZE];

typedef struct IW_register_event_struct *IW_REGISTER_EVENT_PTR;
typedef struct IW_register_event_struct *IW_REGISTER_DIS_CHG_PTR;
typedef struct IW_register_event_struct {
    IVEPtrRep prev_off;
    IVEPtrRep next_off;
    int wid;
    int wave;
    Window send_event_window;
    long event_mask;
    int buttons;
} IW_REGISTER_EVENT,IW_REGISTER_DIS_CHG;

typedef struct IW_WL_ {
    int     prev;         /* previous window struct */
    int     next;         /* next window struct */
    int     scr_x_pos;
    int     scr_y_pos;    /* x and y position in pixels for the upper-left
			     corner of the monitor window (this is the window
			     uppermost in the tree and not necessarily the
			     shell because of window manager reparenting).
			     They are measured relative to the upper-left
			     corner of the screen's root window (positive
			     to the right and down). */
    int     max_width;
    int     max_height;   /* Width and height in pixels of the monitor
			     window (topmost in the tree and not necessarily
			     the shell because of window manager reparenting.
			  */
    int     win_x_pos;
    int     win_y_pos;    /* x and y position in pixels for the drawing
			     widgets' upper-left corner measured relative to
			     the upper-left corner of the window whose
			     size and position is given by max_width,
			     max_height, scr_x_pos, scr_y_pos. */
    int     win_width;
    int     win_height;   /* Width and height in pixels for the drawing
			     widgets. */
    int     modified;     /* Sets whether or not an IW_BANK_DISPLAY_CHANGED
                             event is broadcast on the next redraw (the
                             monitor may temporarily override this while
                             handling repeated mouse/button/slider events
                             and send just one event at the end of the
                             sequence).  One bit (IW_SEC_MODIFIED_BIT) is set
                             when the section number displayed may have
                             changed - i.e. IWAlDisSec or an equivalent
                             action; another bit (IW_MAP_MODIFIED_BIT) is set
                             if there was a possible change in the wave
                             mapping; and another bit (IW_DISP_MODIFIED_BIT)
                             flags any one of a number of other changes:
                             new data written to the section or IWAlZoom,
                             IWAlDisOffset, IWAlMulDisp, IWAlScl, IWAlDisImg
                             and in the appropriate circumstances
                             IWAlSclAlgorithm, IWAlComplexDis, and IWAlStepAttr
                             or any of their equivalents. */
    int     dis_mode;     /* Display mode.   0 : Pseudo Color 1 : True Color */
    int     exist;        /* 0 for nonexist, 1 for exist */
    IVEPtrRep bank_off[IW_MAX_WAVE];
                          /* Offsets to IW_BAL structures holding the bank
			     attributes */
    int     bank_color[IW_MAX_WAVE];
                          /* IW_BLACK, IW_WHITE, IW_RED, IW_GREEN,
			     IW_REDGREEN, IW_REDBLUE, or
			     IW_GREENBLUE */
    IVEPtrRep data_volume_gil_off;
                          /* Offset to the MRC header + extra data */
    IVEPtrRep il_array;   /* Offset to first element of an array of
			   * offsets to IW_IAL structures */
    int     pseudo_graphics_color;
    unsigned long graphics_id;
    IVEPtrRep graphics_list;
                          /* Offset to a list of IW_GRL structures */
    int     graphics_list_size;
    int     graphics_displayed;
    int     gr_z_range[2];
    int     gr_t_range[2];
    int     clear_bkg;
    int     border_style; /* Controls both window manager decorations
			     and application menu/tools shown see
			     IWMacro.h for macros to extract and set the
			     two components. */
    int     md_x_off;
    int     md_y_off;
    int     md_dis_off;

    int     iconified;
    int     delay;
    int     monitor_pid;
    int     monitor_wait;
    XID     monitor_font;
    Window  xwindow[3];
    IVEShmSemaphore sync_enter_lock;
    IVEShmMutex data_sync_lock;
    int     scratch_mode;
    int     cur_wave_num;
    int     mult_disp_x;
    int     mult_disp_y;
    IVEPtrRep scale_gr_off;   /* Offset to IW_GRL describing the scale bar */
    int     scale_bar_vertical;
    float   scale_bar_length;
    float   zoom;         /* Current zoom factor. */
    IVEPtrRep first_reg_event;    /* Offset to IW_REGISTER_EVENT struct */
    IVEPtrRep first_reg_dis_chg;  /* Offset to IW_REGISTER_DIST_CHG struct */
    int     mouse_button[3];
    int     step_type;    
    int     interpolation;
    int     step_inc;
    int     apply_to_wave;
    int     proc_id;
    int     complex_disp;
    int     scaling_algorithm;
    int     img_displayed;
    int     disp_size;
    int     mem_size;
    int     view_flag;
    char    file_name[2 * IW_FILE_NAME_SIZE];      
    int     local_disk;
    int     cur_sec_size;
    int     save_buf_stat;
    /*
     * Holds the offset to a zero-terminated array of offsets to the headers
     * attached to the window.
     */
    IVEPtrRep gr_obj_header_list;
    /*
     * Would combine mouse_button_ext with mouse_button, but clients statically
     * linked against previous versions then would not work with the modified
     * shared memory structure.
     */
    int     mouse_button_ext[2];
    int     initial_delay;
    int     mat_layout;    /* DATA_FMT or TEXT_FMT */
    int     stereo_offset;
    /* Holds the zero-based index of the currently displayed resolution. */
    int     displayed_resolution;
    /* If nonzero, the section indices are displayed on the image. */
    int     section_number_style;
    /*
     * Holds the zero-based index of the color for the section number display
     * (valid values range up to the number of graphics colors minus one.
     */
    int     section_number_color;
} IW_WL;

/*
   These are bits that can be toggled in the "modified" field of the
   IW_WL structure.
*/
#define IW_MODIFIED_CLEAR        (0)
#define IW_SEC_MODIFIED_BIT      (1 << 0)
#define IW_MAP_MODIFIED_BIT      (1 << 1)
#define IW_DISP_MODIFIED_BIT     (1 << 2)

#define IW_ALL_FIELDS               -1
#define IW_BALXoffset               (1 << 4)
#define IW_BALYoffset               (1 << 5)
#define IW_BALImageListPtr          (1 << 6)
#define IW_BALNumImage		    (1 << 7)
#define IW_BALImageListSize         (1 << 8)
#define IW_BALZoom                  (1 << 10)
#define IW_BALMapped                (1 << 11)
#define IW_BALSelected              (1 << 12)
#define IW_BALInterpolation         (1 << 15)
#define IW_BALMovieStart            (1 << 16)
#define IW_BALMovieEnd              (1 << 17)
#define IW_BALMovieInc              (1 << 18)
#define IW_BALMovieDir              (1 << 19)

#define IW_NOT_CONTIGUOUS           (1 << 0)
#define IW_NOT_SAME_DELTA           (1 << 1)
#define IW_NOT_SAME_SIZE            (1 << 2)
#define IW_NOT_SAME_MODE            (1 << 3)
#define IW_NOT_SAME_SCALE           (1 << 4)

/* Points for graphics primitives. */
#define IW_GR_WC_PT                     0x00000001
#define IW_GR_2D_PT                     0x00000002
#define IW_GR_3D_PT                     0x00000004
#define IW_GR_4D_PT                     0x00000008
#define IW_GR_5D_PT                     0x00000010

/* slice_mode in IW_display_image_info */
#define IW_SLICE_XY                      0
#define IW_SLICE_YZ                      1
#define IW_SLICE_XZ                      2
#define IW_SLICE_ANGLE                   3

typedef struct IW_BAL_ {
  int       num_imgs;     /* Number of images currently in IL. */
  int       mapped;       /* BOOLEAN. */
  int       gr_color;
  int       cur_z_num;
  int       cur_time_num;
  int       x_offset;
  int       y_offset;
  int       off_group;
} IW_BAL;

typedef union IW_utag {
  unsigned char byte_value;
  short short_value;
  float float_value;
  IW_COMPLEX_VALUE complex_value;
} IW_UNION_VALUE;

#define IW_IL_IMG_RM_IMG                    0x00000001
#define IW_IL_IMG_RM_BUF                    0x00000010

#define IW_IL_IMG_FREE                      0
#define IW_IL_IMG_MOD_DATA                  1
#define IW_IL_IMG_IMG_DATA                  2
#define IW_IL_IMG_NUM_QUE                   3

#define IW_IL_IMG_IMG_MEM                   IW_IL_IMG_NUM_QUE + 0
#define IW_IL_IMG_DISP_BUF                  IW_IL_IMG_NUM_QUE + 1
#define IW_IL_IMG_NO_COPY                   IW_IL_IMG_NUM_QUE + 2
#define IW_IL_IMG_NOT_LOAD                  IW_IL_IMG_NUM_QUE + 3

#define IW_IL_IMG_EMPTY                    -1
#define IW_IL_IMG_IN_MEM                    0
#define IW_IL_IMG_IN_TMP                    1
#define IW_IL_IMG_IN_FILE                   2
#define IW_IL_IMG_IN_MEM_AND_FILE           3

typedef struct IW_IAL_ {
    float        scale_min;
    float        scale_max;
    float        scale_coeff1;
    float        scale_coeff2;
    float        min;
    float        max;
    float        mean;
    int          buf_stat;
    IVEPtrRep    img_mem_off;
    IVEPtrRep    buf_mem_off;
    IVEPtrRep    img_mem_node;       /* Offset to memStruct for data */
    IVEPtrRep    buf_mem_node;       /* Offset to memStruct for scaled image */
    int          sec_loc;
    IVEPtrRep    graphics_list;      /* Offset to the first IW_GRL structure */
    float        inten_thresh;
    int          inten_thresh_on;
} IW_IAL;

#define MRC_DVID                       -16224

#define MRC_FORMAT                     1
#define MRC_NON_INTERLEAVED            IW_ZTW  
#define MRC_INTERLEAVED                IW_WZT
#define MRC_MAX_WAVES                  5
#define MRC_MULTIPLE_TIMES             6
#define MRC_MULTIPLE_WAVES             7
#define MRC_MULTIPLE_TIMES_AND_WAVES   8

/* Extended header for multiple wavelength data. */
typedef struct IW_MRC_Ext_Header {
  struct {
    float amin;
    float amax;
    float amean;
  } stats[MRC_MAX_WAVES];
} IW_MRC_EXT_HEADER, *IW_MRC_EXT_HEADER_PTR;

/* This is the data volume specification structure. It is for images within
   the same bank with identical geometry specifications. */
typedef struct IW_GIL_ {
  int       x_size;        /*  1       nx */
  int       y_size;        /*  2       ny */
  int       z_size;        /*  3       nz */
  int       mode;          /*  4       mode, Data format in memory. */
  int       ld_x_st;       /*  5       nxstart */
  int       ld_y_st;       /*  6       nystart */
  int       ld_z_st;       /*  7       nzstart */
  int       file_mxyz[3];  /*  8 9 10  mx, my, mz, Default to 1 1 1 */
  float     xyz_len[3];    /* 11 12 13 Pixel spacing in real-world coord. */
/* old def for cel
  float     x_len;         
  float     y_len;
  float     z_len;
*/
  float     file_angle[3]; /* 14 15 16 alpha, beta, gamma to 90.0 90.0 90.0 */
  int       file_mapc;     /* 17       mapc */
  int       file_mapr;     /* 18       mapr */
  int       file_maps;     /* 19       maps */
  float     pad_amin; /* 20 35 37 39 44 */
  float     pad_amax; /* 21 36 38 40 45 */
  float     pad_amean;/* 22 */
  int       file_space_group;       /* 23 ispg: default 0 */
  int       file_ext_header_size;   /* 24 isymb In bytes. */
  union {
      struct {
	  short     nDVID;         /* 25A */
	  short     nblank;        /* 25B */
      } p;
      struct {
	  char extra0[4];
      } e;
  } va;
  int       ld_t_st;       /* 26  */
  char      ibyte[24];     /* 27 - 32 */
  short     file_num_ints;          /* 33A nint */
  short     file_num_reals;         /* 33B nreal */
  short     num_resolutions;        /* 34A nres */
  short     file_nzfact; /* 34B 1:Bining done in x&y. 2:Bining done in x,y,&z.*/
  union {
      struct {
	  float     pad_amin1;         /* 35 */
	  float     pad_amax1;         /* 36 */
	  float     pad_amin2;         /* 37 */
	  float     pad_amax2;         /* 38 */
	  float     pad_amin3;         /* 39 */
	  float     pad_amax3;         /* 40 */
      } p;
      struct {
	  char extra2[24];
      } e;
  } vb;
  short     file_type;     /* 41A file_type */
  short     file_lens;     /* 41B lens default = 8 */
  short     tilt_set_axis; /* 42A n1 for type 1,2 == ND1 */
  short     tilt_set_angle;/* 42B n2 for type 1,2 == VD1 */ 
  short     sects_in_stereo;        /* 43A v1 for type 3,4 == ND1 */
  short     sects_between_pair;     /* 43A v2 for type 3,4 == ND2 */
  union {
      struct {
	  float     pad_amin4;         /* 44 */
	  float     pad_amax4;         /* 45 */
	  short     time_size;         /* 46A time_size for type 6,8 == ND1 */
	  short     wave_focus_order;  /* 46B 0 for noninterleave ; 1 for interleave */
      } p;
      struct {
	  char extra3[12];
      } e;
  } vc;
  float     x_tilt_angle;  /* 47 titl_x Tilt angle from data collection. */
  float     y_tilt_angle;  /* 48 titl_y Tilt angle from data collection. */
  float     z_tilt_angle;  /* 49 titl_z Tilt angle from data collection. */
  union {
      struct {
	  short     wave_size;     /* 50A num_waves */
	  short     wave_length[IW_MAX_WAVE]; /* 50B 51 52 */
	  float     z_orig;        /* 53 Same as above. */
	  float     x_orig;        /* 54 Off from orig. data collect in realworld coor*/
	  float     y_orig;        /* 55 Same as above. */
      } p;
      struct {
	  float x_orig;
	  float y_orig;
	  float z_orig;
	  char map[4];
	  char machst[4];
	  float rms;
      } e;
  } vd;
  int       num_labels;    /* 56 */
  char      label[10][80];
  int       file_dec_format;   /* Boolean. */
  int       nocopy;
  int       hdr_format;
  IVEPtrRep file_ext_header;  /* Offset to extended header in shared memory */
  char      file_name[IW_MAX_WAVE][IW_FILE_NAME_SIZE];      

  /* Data loading params */
  int       valid_ld_info;
  int       ld_x_end, ld_x_bin;
  int       ld_y_end, ld_y_bin;
  int       ld_z_end, ld_z_bin;
  int       ld_z_inc;
  int       ld_time_st, ld_time_end, ld_time_bin;
  int       ld_time_inc;
  int       ld_wave[IW_MAX_WAVE];
  int       scaled;       /* Are images in this DV scaled? */
  float     scale_min[IW_MAX_WAVE];
  float     scale_max[IW_MAX_WAVE];
  float     scale_coeff1[IW_MAX_WAVE];
  float     scale_coeff2[IW_MAX_WAVE];
  int       ld_st_end_erase;
  int       ld_wave_wid[IW_MAX_WAVE];
  int       ld_wave_bank[IW_MAX_WAVE];
  int       ld_wave_bank_pos[IW_MAX_WAVE];

  IW_GIL_PTR prev; /* Pointer to previous element in list. */
  IW_GIL_PTR next; /* Pointer to next element in list. */
  int       num_links;    /* Number of IALs refering to this GIL. */
  char      local_file[IW_FILE_NAME_SIZE];      
  int       set_local_file;
  float     inten_thresh[IW_MAX_WAVE];
  int       inten_thresh_on[IW_MAX_WAVE];
} IW_GIL;

typedef struct IW_shared_block_struct *IW_SHARED_BLOCK_PTR;
typedef struct IW_shared_block_struct {
    IVEPtrRep prev_off;
    IVEPtrRep next_off;
    IVEPtrRep name_off;
    IVEPtrRep data_off;
} IW_SHARED_BLOCK;

typedef struct IW_shared_mem_struct {
    int          compat_major;  /* Only allow a statically-linked application
				   to attach if compat_major matches
				   IW_COMPAT_MAJOR. */
    int          compat_minor;  /* Only allow a statically-linked application
				   to attach if compat_minor is greater than
				   or equal to IW_COMPAT_MINOR. */
    IVEPtrRep    WL_array;
    int          max_window; /* max number of windows that are allocated */
    int          top_window; /* top window in the list. */
    int          bot_window; /* bot window in the list. */
    int          num_window; /* number of windows created */ 
    IVEShmMutex  focus_lock;
    int          focus_window;
    int          num_displayed; /* number of windows displayed. */
    int          screen_left;
    int          screen_bottom;
    int          screen_width;
    int          screen_height;
    int          scale_max;   /* max color index used to display data */
    int          scale_min;   /* minimum color index used to display data */
    int          graphics_system;

    Colormap     pseudo_colormap;

    Window       ind_bk_win;

    IVEPtrRep    top_shared_block;    /* Offset of the first IW_SHARED_BLOCK */

    int          graphics_id;
    int          rgb_overlay_color[IW_NUM_COLORS];
    long         img_mem_size;
    long         total_img_size;
    long         mem_mode;
    IVEShmHold   mem_lock;
    memQueue     mem_queue[IW_IL_IMG_NUM_QUE];
    int          que_pri[IW_IL_IMG_NUM_QUE];
    long         pg_thresh;
    long         num_obj_lists;
    IVEPtrRep    obj_lists;
    IVEShmMutex  block_lock;  /* Lock for top_shared_block and children. */
} IW_SHARED_MEM_STRUCT;

typedef struct IW_stuff_ {
  int ispg;
  int inbsym;
  char room[32];
  short nint;
  short nreal;
  short nres;
  short nzfact;
  float  wmaxmin[6];
  short data_type;
  short lens;
  short n1;
  short n2;
  short v1;
  short v2;
  float w5maxmin[2];
  short ntime;
  short wzt;
  float tilt[3];
  short nwaves;
  short wavelength[5];
  float zorig; 
} IW_STUFF;

typedef struct IW_imsubs_ {
  int lstream[24];
  int nbhdr;
  int nbw;
  int nbw3;
  int nb[5];
  int nbl;
  int ncr3[10][3];
  int mode[10];
  int ncrst[10][3];
  int nxyz[10][3];
  float cel[10][6];
  int mapcrs[10][3];
  float denmmm[10][3];
  IW_STUFF stuff[10];
  float origxy[10][2];
  int nlab[10];
  char label[10][800];
  int flag[10];
  int nocon[10];
  int nbsym[10];
  int ibsym[10][15000];
  int npointer[10][5];
  int numopen;
  int print;
  int window[10];
} IW_IMSUB, *IW_IMSUB_PTR;

typedef struct window_bank_element {
  int window;
  int bank;
} IW_WIN_BANK, *IW_WIN_BANK_PTR;

typedef IW_IAL                       IW_IMAGE_INFO;
typedef IW_IAL_PTR                   IW_IMAGE_INFO_PTR;
typedef IW_GIL                       IW_IMAGES_INFO;
typedef IW_GIL_PTR                   IW_IMAGES_INFO_PTR;
typedef int                          BOOLEAN;

/*
 * Used to pass control information from the client to the window for
 * IWCaptureImage.
 */
typedef struct IWCaptureControl {
    IWEncodedSHMPtr buffer_offset;
    int ix, iy, width, height;
    int buffers, format, type;
    int status;
} IWCaptureControl;

#endif /* include guard */
