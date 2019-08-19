#ifndef UCSF_MSG_LUT_H
#define UCSF_MSG_LUT_H

/*
 * Declares constants, structures, and function prototypes for manipulating
 * lookup tables for image colors.
 */

#define LUT_FMT_INVALID (1)
#define LUT_FMT_PAL (2)
#define LUT_FMT_VTK (3)

#define LUT_OPT_REVERSE (1)
#define LUT_OPT_SCALE (2)
#define LUT_OPT_SHIFT (4)
#define LUT_OPT_SKIP (8)
#define LUT_OPT_TRUNCATE (16)

#define LUT_ERR_GENERIC (1)
#define LUT_ERR_MEM (2)
#define LUT_ERR_FILE (3)
#define LUT_ERR_IO (4)
#define LUT_ERR_SYNTAX (5)
#define LUT_ERR_OVERFLOW (6)
#define LUT_ERR_ARG (7)
#define LUT_ERR_UNSUPPORTED (8)
#define LUT_ERR_INTERNAL (9)

struct LUTListImpl;
typedef struct LUTListImpl* LUTList;

struct LUTProxyImpl;
typedef struct LUTProxyImpl* LUTProxy;

typedef struct LUTFileAttr {
    int format;
    int fixed;
    union {
	struct { int minsize, maxsize; } variable;
	struct { int size; } fixed;
    } v;
} LUTFileAttr;

#ifdef __cplusplus
extern "C" {
#endif

/* Creates an empty list of lookup tables. */
LUTList lut_create_list(void);

void lut_destroy_list(LUTList list);

/*
 * Appends the lookup tables listed in a file to an existing list.
 * Returns zero if successful; otherwise returns a nonzero value and
 * leaves the input list intact.  The value of desired_size acts as
 * a filter on the color tables extracted from the file.  If desired_size
 * is positive, only those tables that can match that size are extracted.  If
 * desired_size is zero, all tables are extracted.
 */
int lut_append_to_list_from_file(
    const char* filename, int desired_size, int error_ignored, LUTList list
);

/* Returns zero if successful and a nonzero value if not successful. */
int lut_save_list(LUTList list, const char* filename);

int lut_get_list_count(LUTList list);

LUTProxy lut_get_list_element(LUTList list, int index);

/*
 * Inserts the given lookup table before the position index in the
 * list.  Returns zero if successful and a nonzero value if not successful.
 * The list assumes responsibility for freeing any resources associated with
 * the inserted element.
 */
int lut_insert_list_element(LUTList list, int index, LUTProxy proxy);

void lut_delete_list_range(LUTList list, int start, int count);

LUTProxy lut_create_proxy(
    const char* name, const char* short_name, const char* path, int options
);

void lut_destroy_proxy(LUTProxy proxy);

const char* lut_get_proxy_name(LUTProxy proxy);

const char* lut_get_proxy_short_name(LUTProxy proxy);

const char* lut_get_proxy_path(LUTProxy proxy);

int lut_get_proxy_options(LUTProxy proxy);

int lut_set_proxy_name(LUTProxy proxy, const char* name);

int lut_set_proxy_short_name(LUTProxy proxy, const char* short_name);

int lut_set_proxy_path(LUTProxy proxy, const char* path);

int lut_set_proxy_options(LUTProxy proxy, int options);

int lut_lookup_rgb_int_from_proxy(
    LUTProxy proxy, int n, int index, int* p_red, int* p_green, int* p_blue
);

int lut_check_file_attr(const char* path, LUTFileAttr* p_attr);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
