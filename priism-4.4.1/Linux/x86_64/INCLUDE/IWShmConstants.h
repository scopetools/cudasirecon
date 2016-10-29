#ifndef UCSF_MSG_IWSHMCONSTANTS_H
#define UCSF_MSG_IWSHMCONSTANTS_H

#include <sys/types.h>  /* pid_t */
#include "ive_shm.h"    /* IVEPtrRep */

typedef struct IW_GRL_ *IW_GRL_PTR; /* Graphics List. */

typedef struct {
    float x;
    float y;
} IW_POINT_2D, *IW_POINT_2D_PTR;

typedef struct {
    float x;
    float y;
    float z;
} IW_POINT, *IW_POINT_PTR;

typedef struct IW_3D_Point{
    float   coord[3];
    float   normal[3];
    float   radius;
    int     marker;
    int     pol_cnt;
    IVEPtrRep label_off;
    IVEPtrRep prev_off;     /* offset to a IW_3D_Point */
    IVEPtrRep next_off;     /* offset to a IW_3D_Point */
    IVEPtrRep prev_sib_off; /* offset to an IW_3D_Point */
    IVEPtrRep next_sib_off; /* offset to an IW_3D_Point */
} IW_3D_Point;

typedef struct IW_GRL_ {
    IVEPtrRep prev_off;
    IVEPtrRep next_off;
    int displayed;
    int marked_for_delete;
    int primitive_type;
    IVEPtrRep primitive_off;
} IW_GRL;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
} IW_GR_ANY, *IW_GR_ANY_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    float  x;
    float  y;
    float  z;
    int    shape;
} IW_GR_POINT, *IW_GR_POINT_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    float  x1;
    float  y1;
    float  z1;
    float  x2;
    float  y2;
    float  z2;
    int    thickness;
    int    pattern;
} IW_GR_LINE, *IW_GR_LINE_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    float  x;
    float  y;
    float  z;
    float  width;
    float  height;
    float  depth;
    int    thickness;
    int    pattern;
} IW_GR_BOX, *IW_GR_BOX_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    float  x;
    float  y;
    float  z;
    float  radius;
    int    thickness;
} IW_GR_CIRCLE, *IW_GR_CIRCLE_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    float  x;
    float  y;
    float  z;
    int    size;
    IVEPtrRep string_off;
} IW_GR_STRING, *IW_GR_STRING_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    IVEPtrRep points_list;  /* Offset to first IW_POINT */
    IVEPtrRep color_list;   /* Offset to first color (int) */
    IVEPtrRep string_list;  /* Offset to first string (all stored contig.) */
    int    max_char;
    int    size;
    int    num_strings;
} IW_GR_MULT_STRING, *IW_GR_MULT_STRING_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    IVEPtrRep points_list;   /* Offset to first IW_POINT */
    int    num_points;
    int    shape;
} IW_GR_POINTS, *IW_GR_POINTS_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    IVEPtrRep points_list;   /* Offset to first IW_POINT */
    int    num_points;
    int    thickness;
    int    pattern;
} IW_GR_LINES, *IW_GR_LINES_PTR;

#define IW_FILL_POLYGON          -1
#define IW_MAX_POLYGON_SIZE      256
typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    IVEPtrRep points_list;   /* Offset to first IW_POINT */
    int    num_points;
    int    thickness;  /* or IW_FILL_POLYGON */
    int    pattern;
} IW_GR_POLYGON, *IW_GR_POLYGON_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    IVEPtrRep points_list;   /* Offset to first IW_POINT */
    int    num_points;
    int    thickness;
    int    pattern;
} IW_GR_MULT_LINES, *IW_GR_MULT_LINES_PTR;

typedef struct {
    int    id;
    int    color;
    int    displayed;
    int    type;
    pid_t  gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    float  x1;
    float  y1;
    float  z1;
    float  x2;
    float  y2;
    float  z2;
    int    thickness;
    int    pattern;
} IW_GR_3D_LINE, *IW_GR_3D_LINE_PTR;

typedef struct JustaDot{
    float   coord[3];
    float   normal[3];
    IVEPtrRep nextdot_off;  /* offset to JustaDot */
}JustaDot;

typedef struct IW_3D_List{
    struct  IW_3D_Point p_list[50];
    IVEPtrRep next_list;    /* offset to an IW_3D_list */
} IW_3D_List;

typedef struct IW_GR_3D_MODEL{
    int    id;
    int    color;
    int    displayed;
    int    type;
    int    gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    int    obj_num;
    int    thickness;
    int    number_of_points;
    float  def_radius;
    float  check_pt[3];
    IVEPtrRep name_off;
    IVEPtrRep first_off;          /* offset to an IW_3D_Point */
    IVEPtrRep last_off;           /* offset to an IW_3D_Point */
    IVEPtrRep current_point_off;  /* offset to an IW_3D_Point */
    struct IW_3D_List point_list;
    IVEPtrRep dotlist;            /* offset to a JustaDot */
    IVEPtrRep lastdot;            /* offset to a JustaDot */
    IVEPtrRep prev_off;           /* offset to an IW_GR_3D_MODEL */
    IVEPtrRep next_off;           /* offset to an IW_GR_3D_MODEL */
} IW_GR_3D_MODEL;

typedef struct IW_GR_5D_MODEL{
    int    id;
    int    color;
    int    displayed;
    int    type;
    int    gr_id;
    int    time_num;
    int    wave_num;
    int    dummy3;
    int    dummy2;
    int    dummy1;
    int    last_obj;
    IVEPtrRep current_obj_off;  /* offset to an IW_GR_3D_MODEL */
    IVEPtrRep bottom_obj_off;   /* offset to an IW_GR_3D_MODEL */
    IVEPtrRep top_obj_off;      /* offset to an IW_GR_3D_MODEL */
    IVEPtrRep prev_off;         /* offset to the previous IW_GR_5D_MODEL */
    IVEPtrRep next_off;         /* offset to the next IW_GR_5D_MODEL */
} IW_GR_5D_MODEL, *IW_GR_5D_MODEL_PTR;


typedef struct {
    float act_point[3];
    IVEPtrRep mod_list;         /* offset to IW_GR_5D_MODEL */
    int cur_time;
    int display_all;
} IW_MODEL_UPDATE_INFO;

#define IW_GR_Point                     0x00000001
#define IW_GR_Line                      0x00000002
#define IW_GR_Box                       0x00000004
#define IW_GR_Circle                    0x00000008
#define IW_GR_String                    0x00000010
#define IW_GR_Points                    0x00000020
#define IW_GR_Lines                     0x00000040
#define IW_GR_Polygon                   0x00000080
#define IW_GR_Mult_Lines                0x00000100
#define IW_GR_3D_Line                   0x00000200
#define IW_GR_3D_Model                  0x00000400
#define IW_GR_Mult_Strings              0x00000800

#define IW_GR_DISPLAYED			0x00000001
#define IW_3D_DISP_ALWAYS		0x00000002
#define IW_3D_DISP_NORMAL		0xFFFFFFFD
#define IW_3D_DETACH                    0x00000004
#define IW_GR_DEL_MARKED		0x00000008
#define IW_GR_DEL_ALL   		0x00000010
#define IW_GR_DEL_PROC  		0x00000020
#define IW_GR_DEL_RESET 		0xFFFFFFE7

/* Patterns for graphics primitives. */
#define IW_GR_SOLID                      0
#define IW_GR_DASHED                     1
#define IW_GR_MORSE                      2

/* Points for graphics primitives. */
#define IW_GR_DOT                        0
#define IW_GR_X                          1
#define IW_GR_O                          2
#define IW_GR_PLUS                       3
#define IW_GR_SQUARE                     4
#define IW_GR_TRIANGLE                   5
#define IW_GR_DIAMOND                    6

/*
 * Create an alias for IVEPtrRep for users of IWEncodeSHMPtr and
 * IWDecodeSHMPtr.
 */
typedef IVEPtrRep IWEncodedSHMPtr;

#endif /* include guard */
