#ifndef UCSF_MSG_PICKFUNCS_H
#define UCSF_MSG_PICKFUNCS_H

/*
 * A Picklist is a list of points in a 3D volume.  Each point has four
 * predefined attributes: x, y, and z coordinates and an integer "name".
 * Additional attributes can be added to a point dynamically.  Attributes
 * are also assignable to the list.
 *
 * A Pickset is a group of Picklists referenced by wavelength number,
 * time point number, and list number.  Locking and I/O operations are
 * available at the granularity of a Pickset.
 */

#include "shardict.h"
#include "privdict.h"
#include <limits.h>          /* INT_MIN */
#include <stdio.h>           /* FILE */


#define PICKSET_INVALID_INDEX (INT_MIN)


struct PicklistImpl;
typedef struct PicklistImpl* Picklist;
typedef const struct PicklistImpl* PicklistRO;

typedef struct PicklistIter {
    Picklist list;
    int index;
} PicklistIter;
typedef struct PicklistIterRO {
    PicklistRO list;
    int index;
} PicklistIterRO;

struct PicksetImpl;
typedef struct PicksetImpl* Pickset;
typedef const struct PicksetImpl* PicksetRO;

typedef struct PicksetIter {
    union {
	PrivdictIter p;
	ShardictIter s;
    } var;
    int shared;
} PicksetIter;
typedef struct PicksetIterRO {
    union {
	PrivdictIterRO p;
	ShardictIterRO s;
    } var;
    int shared;
} PicksetIterRO;

typedef struct Pickpt {
    float x, y, z;
    int name;
} Pickpt;

struct PicksetSaveCustomizationImpl;
typedef struct PicksetSaveCustomizationImpl* PicksetSaveCustomization;

struct PicksetNTWOrderImpl;
typedef struct PicksetNTWOrderImpl* PicksetNTWOrder;

typedef int (*PicksetSaveCustomHeader)(PicksetRO pickset, void* arg, FILE* fo);
typedef int (*PicksetSaveInitIterator)(PicksetRO pickset, void* arg);
typedef int (*PicksetSaveIncrIterator)(
    PicksetRO pickset, void* arg, FILE* fo, PicksetIterRO* p_iter
);
typedef void (*PicksetSaveCleanIterator)(PicksetRO pickset, void* arg);
typedef int (*PicksetSaveCustomListHeader)(
    PicksetRO pickset, void* arg, FILE* fo, PicksetIterRO list_iter
);
typedef int (*PicksetSaveCustomPoint)(
    PicksetRO pickset,
    void* arg,
    FILE* fo,
    PicksetIterRO list_iter,
    PicklistIterRO pt_iter
);

typedef int (*PicksetPrintWTInfo)(
    PicksetRO pickset, void* arg, FILE* fo, int iwave, int itime
);

struct PicksetLoadCustomizationImpl;
typedef struct PicksetLoadCustomizationImpl* PicksetLoadCustomization;

typedef int (*PicksetLoadCustomComment)(
    int line_number, char* line, Pickset pickset, void* arg, int* p_offset
);
typedef int (*PicksetLoadCustomWave)(
    int line_number,
    char* line,
    Pickset pickset,
    void* arg,
    int* p_offset,
    int* p_wave,
    int* p_time
);
typedef int (*PicksetLoadCustomList)(
    int line_number,
    char* line,
    Pickset pickset,
    void* arg,
    int* p_offset,
    int iwave,
    int itime,
    PicksetIter* p_iter
);
typedef int (*PicksetLoadCustomPoint)(
    int line_number,
    char* line,
    Pickset pickset,
    void* arg,
    int* p_offset,
    PicksetIter list_iter,
    PicklistIter* p_iter
);
typedef int (*PicksetLoadCustomMisc)(
    int line_number, char* line, Pickset pickset, void* arg, int* p_offset
);

/*
 * im_stream_num is the IM library stream number which is attached to an
 * open image stack.  filename is a pointer to a null-terminated character
 * string which is the name of the image stack.  filename may be NULL.
 */
typedef struct PicksetSourceData {
    int im_stream_num; char* filename;
} PicksetSourceData;

/*
 * flags is a bitwise or of 0, PICKSET_WAVE_UNINIT_ZERO,
 * PICKSET_TIME_UNINIT_ZERO, PICKSET_WAVE_ALLOW_OBOUND,
 * PICKSET_TIME_ALLOW_OBOUND.
 */
enum {
    PICKSET_WAVE_UNINIT = 1,
    PICKSET_TIME_UNINIT = 2,
    PICKSET_WAVE_OBOUND = 4,
    PICKSET_TIME_OBOUND = 8
};
typedef struct PicksetWTBounds {
    const int* waves_used;
    int wave_count;
    int wave_offset;
    int wave_uninit;
    int time_min;
    int time_max;
    int time_step;
    int time_uninit;
    int flags;
} PicksetWTBounds;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Creates a new Pickset which is private (visible only to the application).
 *
 * Returns a non-NULL value if the function call succeeded and a NULL value
 * When no longger needed, the resources associated with the Pickset can
 * be freed with pickset_detach or pickset_delete.
 */
Pickset pickset_create(void);

/*
 * Creates a new shared (visible to this application and others) Pickset
 * or attaches to an existing one.  name is a null-terminated string with
 * the name of the set.  There is a standard name used for the Pickset
 * associated with a particular image window; pickset_fill_name can be used
 * to generate such a name.  flags must be a bitwise or of one or more of the
 * following: 0, PICKSET_CREATE, PICKSET_EXCL, PICKSET_CLEAR, and
 * PICKSET_WLOCK.  If PICKSET_CREATE is set, a new Pickset will be created
 * if necessary.  If PICKSET_EXCL is set, the function call will fail if an
 * Pickset with the same name already exists.  If PICKSET_CLEAR is set and the
 * call attaches to an existing Pickset, the Pickset will be cleared before
 * pickset_attach_shared returns.  If PICKSET_WLOCK is set, the list will
 * be locked for writing when pickset_attach_shared returns; it is the
 * caller's responsibility to release the lock when it is no longer needed.
 *
 * Returns a non-NULL value if the function call succeeded and a NULL
 * value if it failed.  If the application wants to determine the cause
 * of failure, it should pass a non-NULL value for p_status.  Then, if
 * the call fails, *p_status will be set to one of the following:
 *
 * IVE_ERR_NO_MEM
 *      There was insufficient space or system resources.
 * IVE_ERR_EXIST
 *      PICKSET_EXCL was set but an array with the same name was found.
 * IVE_ERR_INVALID
 *      PICKSET_CREATE was not set and there was no array to attach to,
 *
 * pickset_attach_shared may deadlock if the caller holds a read or write
 * lock to the array with the given name.
 */
enum {
    PICKSET_CREATE = 1, PICKSET_EXCL = 2, PICKSET_CLEAR = 4, PICKSET_WLOCK = 8
};
Pickset pickset_attach_shared(int flags, const char* name, int* p_status);

/*
 * Given a character array, name, that must be at least PICKSET_NAME_SIZE,
 * characters long, fills it with standard name for the Pickset associated
 * with the image window whose number is window_id.
 */
#define PICKSET_NAME_SIZE 30
void pickset_fill_name(int window_id, char* name);

/*
 * Detaches from a Pickset.  if this is the last application attached to
 * the Pickset (which is always the case for Picksets created with
 * pickset_create), the resources associated with the list are released.
 *
 * For a shared list, picklist_detach may not be called if the caller has
 * locked the list.
 */
void pickset_detach(Pickset pickset);

/*
 * Deletes the resources for a Pickset.  If the set of lists is shared, this
 * may be dangerous (there may be other applications still using the list).
 *
 * pickset_delete may not be called if the caller has locked the set of lists.
 */
void pickset_delete(Pickset pickset);

/*
 * Lock or unlock the set of lists for reading or writing.  These calls have no
 * effect on a private set of lists.
 */
void pickset_lock_read(Pickset pickset);
void pickset_lock_write(Pickset pickset);
void pickset_unlock_read(Pickset pickset);
void pickset_unlock_write(Pickset Pickset);

/*
 * Returns true if the set of lists is shared and false if not.
 */
int pickset_shared(PicksetRO pickset);

/*
 * Returns the number of lists in a set of lists.
 */
int pickset_size(PicksetRO pickset);

/*
 * Removes all lists from a Pickset.
 */
void pickset_clear(Pickset pickset);

/*
 * Saves an set of lists to an ASCII file.  Certain aspects of the output
 * may be customized by passing a non-NULL custom parameter.  See
 * pickset_create_custom_save for details and what the default behavior is
 * when a NULL custom parameter is used.
 *
 * Returns IVE_ERR_NONE if successful or one of the following if not (other
 * values are possible if the caller customizes the output and one of the
 * customization functions returns a value that is not IVE_ERR_NONE or one of
 * the following):
 *
 *      IVE_ERR_BAD_FILENAME
 *           file_name could not be opened for writing
 *      IVE_ERR_IO
 *           An I/O error occurred while writing the file.
 *      IVE_ERR_INVALID
 *           pickset was invalid.
 */
int pickset_save(
    PicksetRO pickset,
    const char* file_name,
    PicksetSaveCustomization custom
);

/*
 * Creates an object which can be passed to pickset_save to customize
 * what is done when the set of points lists is saved.  The object returned
 * specifies the default behavior:
 *
 *     No header is printed.
 *
 *     The lists are printed in the default ordering
 *
 *     The line that in the file that introduces the list has nothing
 *     printed after the leading 'L'.
 *
 *     For each point the x, y, and z coordinates are printed in the default
 *     fixed point format but no additional data is printed on the line.
 *
 * To change this make one or more calls to the pickset_add_saveproc_*
 * pickset_set_saveproc_* functions on the customization object.
 *
 * Returns a non-NULL value if success or NULL if there was insufficient
 * memory.  A non-NULL value should be deleted with
 * pickset_delete_custom_save when it is no longer needed.
 */
PicksetSaveCustomization pickset_create_custom_save(void);
void pickset_delete_custom_save(PicksetSaveCustomization custom);

/*
 * Adds a function to the list of functions that is called by pickset_save
 * immediately after the output file is opened.  The functions are called
 * in the order they are registered until there are no more or until a
 * function returns a value not equal to IVE_ERR_NONE.  In case of the latter,
 * the partially written file is deleted and the value is returned to the
 * caller of pickset_save.  The function must accept three arguments:
 *
 *     a PicksetRO which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_saveproc_header is called
 *
 *     a FILE* which is the handle for the output file
 *
 * By convention the lines in the header begin with a '#'.
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_saveproc_header(
    PicksetSaveCustomization custom,
    PicksetSaveCustomHeader func,
    void* arg
);

/*
 * This is a function, compatible with pickset_add_saveproc_header, which
 * assumes that arg is a pointer to a PicksetSourceData structure and prints
 * out information from that IM file's header.  It does not access the
 * Pickset argument and only uses the filename field of the
 * PicksetSourceData structure if the image stack is not an image window.
 *
 * The common use for this function is to provide information at the start
 * a point list file about the data set from which the points were drawn.
 */
int pickset_print_source_data_info(PicksetRO pickset, void* arg, FILE* fo);

/*
 * Sets the functions that are called by pickset_save to determine the next
 * point list to be written to the output file.  init_func is called
 * to initialize any internal state that is needed.  It should return an
 * integer:  IVE_ERR_NONE if successful or some other value which will returned
 * to the caller of pickset_save.  init_func is passed two arguments:
 *
 *     a PicksetRO which is the set of lists being operated upon
 *
 *     a void* which is the value passed as init_arg when
 *         pickset_set_saveproc_order is called
 *
 * incr_func is called for each iteration to determine the list to be saved
 * in that iteration.  iter_func is passed four arguments:
 *
 *     a PicksetRO which is the set of lists being operated upon
 *
 *     a void* which is the value passed as incr_arg when
 *         pickset_set_saveproc_order is called
 *
 *     a FILE* which is the handle for the output file
 *
 *     an PicksetIterRO* which points to an iterator to the last list printed
 *     (or an invalid iterator if this is the first list to be printed); the
 *     function should change the iterator to point to the next list to be
 *     printed
 *
 * incr_func should return an integer value.  pickset_save examines the
 * value returned.  When the return value is IVE_ERR_NONE, pickset_save
 * examines the PicksetIterRO* and prints the corresponding list if the
 * PicksetIterRO* points to a valid iterator or closes the file if the
 * PicksetIterRO* points to an invalid iterator.  If the returned value is
 * something other than IVE_ERR_NONE, the partially written file is deleted
 * and the returned value is passed along to the caller of pickset_save.
 *
 * clean_func is called when iterations are complete to undo the work of
 * init_func.  It does not return a value and is passed two arguments:
 *
 *     a PicksetRO which is the set of lists being operated upon
 *
 *     a void* which is the value passed as clean_arg when
 *         pickset_set_saveproc_order is called
 *
 * Returns zero if successful or one if there was insufficient memory to
 * set the function.
 */
int pickset_set_saveproc_order(
    PicksetSaveCustomization custom,
    PicksetSaveInitIterator init_func,
    void* init_arg,
    PicksetSaveIncrIterator incr_func,
    void* incr_arg,
    PicksetSaveCleanIterator clean_func,
    void* clean_arg
);

/*
 * This is a function, compatible with pickset_set_saveproc_order,
 * which causes the lists to be printed in order of wavelength index (the
 * outermost loop), time point index, and list number (the innermost loop)
 * at the start of each wavelength/time point pair a line of the form
 *
 *     W # T #
 *
 * is printed to the output file.
 *
 * To use, pass pickset_get_next_by_ntw as the incr_func argumet to
 * pickset_set_saveproc_order.  The other arguments to
 * pickset_set_saveproc_order should be:
 *
 * init_func   0
 * init_arg    anything (not used)
 * incr_arg    a derefenceable points to a PicksetNTWOrder
 *             (see pickset_create_ntw_order)
 * clean_func  0
 * clean_arg   anything (not used)
 */
int pickset_get_next_by_ntw(
    PicksetRO pickset, void* arg, FILE* fo, PicksetIterRO* p_iter
);

/*
 * Allocates and intializes a PicksetNTWOrder for use with
 * pickset_get_next_by_ntw and friends.  A time point with index, i, is
 * included in the output if time_min <= i <= time_max and (i - time_min) %
 * time_step == 0.  A wavelength with index, i, is included in the output if
 * 0 <= i < wave_count and waves_used[i] is true.
 *
 * Returns a non-NULL value if successful.  When the returned value is no
 * longer needed, it should be deleted with pickset_delete_ntw_order.
 */
PicksetNTWOrder pickset_create_ntw_order(
    int time_min,
    int time_max,
    int time_step,
    int wave_count,
    const int waves_used[]
);
void pickset_delete_ntw_order(PicksetNTWOrder ntw_order);


/*
 * Adds a function to the list of functions that is called
 * by pickset_get_next_by_ntw to print out the line that introduces a new
 * wavelength/time point pair (pickset_get_next_by_ntw will print the
 * leading W # T # and trailing newline).  The functions are called in the
 * order registered until there are no more or until a function returns a
 * value not equal to IVE_ERR_NONE.  In case of the latter, the value
 * is returned to the caller of pickset_get_next_by_ntw.  The function
 * is passed five arguments:
 *
 *     a PicksetRO which is the set of lists being written
 *
 *     a void* which is the value passed as arg to pickset_create_ntw_order
 *
 *     a FILE* which is the handle to the file being written
 *
 *     an int which is the current wavelength index
 *
 *     an int which is the current time point index
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_saveproc_wt(
    PicksetNTWOrder ntw_order, PicksetPrintWTInfo func, void* arg
);

/*
 * Adds a function to the list of functions that is called by pickset_save to
 * print out the line that introduces a new list (pickset_save will print the
 * leading "L" and trailing new line).  The functions are called in the order
 * registered until there are no more or until a function returns a value not
 * equal to IVE_ERR_NONE.  In case of the latter, the partially written
 * file is deleted and the value is returned to the caller of pickset_save.
 * The function is passed four arguments:
 *
 *     a PicksetRO which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_saveproc_list is called
 *
 *     a FILE* which is the handle for the output file
 *
 *     a PicksetIterRO which is an iterator to the list to be printed
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_saveproc_list(
    PicksetSaveCustomization custom,
    PicksetSaveCustomListHeader func,
    void* arg
);

/*
 * Adds a function to the list of functions that is called by pickset_save to
 * print out the line that with information about a single point
 * (pickset_save will print the trailing new line).  The
 * functions are called in the order they are registered until there are no
 * more functions or a function returns a value not equal to IVE_ERR_NONE.  In
 * case of the latter, the partially written file is deleted and the value is
 * returned to the caller of pickset_save.  The functions is passed five
 * arguments:
 *
 *     a PicksetRO which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_saveproc_list is called
 *
 *     a FILE* which is the handle for the output file
 *
 *     a PicksetIterRO which points to the list current being printed
 *
 *     a PicklistIterRO which points to the point about to be printed
 *
 * The function should return an integer value.  pickset_save examines the
 * value returned.  If it is something other than IVE_ERR_NONE, the partially
 * written file is deleted and the returned value is passed along to the
 * caller of pickset_save.
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_saveproc_point(
    PicksetSaveCustomization custom, PicksetSaveCustomPoint func, void* arg
);

/*
 * This is a function compatible with pickset_add_saveproc_point which
 * prints the coordinates of the point in a %8.3f %8.3f %8.3f format.
 */
int pickset_print_pt_coord(
    PicksetRO pickset,
    void* arg,
    FILE* fo,
    PicksetIterRO list_iter,
    PicklistIterRO pt_iter
);

/*
 * This is a function compatible with pickset_add_saveproc_point which
 * prints the coordinates of the point and the point name in a
 * %8.3f %8.3f %8.3f %d format.
 */
int pickset_print_pt_coord_name(
    PicksetRO pickset,
    void* arg,
    FILE* fo, 
    PicksetIterRO list_iter,
    PicklistIterRO pt_iter
);

/*
 * Reads an set of lists from an ASCII file.  If append is true, the
 * addtional lists will be appended to what is already present, otherwise the
 * current contents are cleared before the new lists are loaded.  Certain
 * aspects of the input may be customized by passing a non-NULL custom
 * parameter.  See pickset_create_custom_load for details and what the
 * default behavior is when a NULL custom parameter is used.  If p_line is not
 * NULL, *p_line will be set to the line number in the file where processing
 * stopped (line 0 indicates before the start of the file).
 *
 * Returns IVE_ERR_NONE if successful and one of the following if not (other
 * values are possible if the caller customizes the output and one of the
 * customization functions returns a value that is not IVE_ERR_NONE or one
 * of the following):
 *
 *      IVE_ERR_BAD_FILENAME
 *           file_name could not be opened for reading.
 *      IVE_ERR_IO
 *           An I/O error occurred while reading the file.
 *      IVE_ERR_CORRUPT
 *           The file did not have the correct format.
 *      IVE_ERR_NO_MEM
 *           There was insufficient space to store the loaded points.
 *      IVE_ERR_INVALID
 *           pickset was invalid.
 *
 */
int pickset_load(
    Pickset pickset,
    const char* file_name,
    int append,
    PicksetLoadCustomization custom,
    int* p_line
);

/*
 * Creates an object which can be passed to pickset_load to customize
 * what is done when the set of points lists is loaded.  The object returned
 * specifies the default behavior:
 *
 *     Any lines beginning with '#' are ignored.
 *
 *     Any line starting with 'W' is assumed to have the form:
 *
 *         W %d T %d[whitespace][anything]
 *
 *     The two integers are scanned for the wavelength and time index to
 *     apply to all lists until the next line starting with 'W'; the rest
 *     of the line is ignored.  A line starting with 'W' which does not
 *     match the format is treated as an error.
 *
 *     A line starting with 'L' introduces a new list.  If there has
 *     been a previous line starting with 'W', the new list is assigned to
 *     that wavelength and time index; otherwise, it is assigned to wavelength
 *     zero and time point zero.  Anything after the leading 'L' is ignored.
 *
 *     Lines following a line starting with 'L' which do not start with a 'W'
 *     or '#' (which are handled as described above) are scanned to read
 *     a point's coordinates.  If in addition, there is a fourth entry on the
 *     the line, a scan is attempted to read it as the integer name of the
 *     point.  Lines which do not start with three values that can be parsed
 *     as floating-point values are treated as an error.
 *
 * To change this make one or more calls to the pickset_add_loadproc_*
 * functions on the customization object.
 *
 * Returns a non-NULL value if success or NULL if there was insufficient
 * memory.  A non-NULL value should be deleted with
 * pickset_delete_custom_load when it is no longer needed.
 */
PicksetLoadCustomization pickset_create_custom_load(void);
void pickset_delete_custom_load(PicksetLoadCustomization custom);

/*
 * Adds a function to the list of functions that is called by pickset_load
 * when a line beginning with '#' is encountered.  The functions are called
 * in the order they are registered until there are no more or *p_offset is
 * set to something less than zero or to the something greater than or equal
 * to the number of characters in the current line.  pickset_load looks at
 * the integer returned by the last function called and if other than
 * IVE_ERR_NONE, aborts the load and returns that value to the caller of
 * pickset_load.  The function must accept five arguments:
 *
 *     an int which is the line number currently being processed
 *
 *     a char* which is the NULL terminated string holding the contents of
 *     the current line
 *
 *     a Pickset which is the array being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_loadproc_comment is called
 *
 *     an int* which points to the offset from the start of the line where
 *     the routine should start processing and should be set before the
 *     function returns to the offset where the next routine should start
 *     processing (or a negative value or value greater than equal to the
 *     length of the line to terminate processing)
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_loadproc_comment(
    PicksetLoadCustomization custom,
    PicksetLoadCustomComment func,
    void* arg
);

/*
 * Adds a function to the list of functions that is called by pickset_load
 * when a line beginning with 'W' is encountered.  The functions are called
 * in the order they are registered until there are no more or *p_offset is
 * set to something less than zero or to the something greater than or equal
 * to the number of characters in the current line.  pickset_load looks at
 * the integer returned by the last function called and if other than
 * IVE_ERR_NONE, aborts the load and returns that value to the caller of
 * pickset_load.  The function must accept seven arguments:
 *
 *     an int which is the line number currently being processed
 *
 *     a char* which is the NULL terminated string holding the contents of
 *     the current line
 *
 *     a Pickset which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_loadproc_wave is called
 *
 *     an int* which points to the offset from the start of the line where
 *     the routine should start processing and should be set before the
 *     function returns to the offset where the next routine should start
 *     processing (or a negative value or value greater than equal to the
 *     length of the line to terminate processing)
 *
 *     an int* which points to the wavelength index; at the start of
 *     processing the line this will be the previous wavelength index if set
 *     previously or PICKSET_INVALID_INDEX if not; may be modified
 *
 *     an int* which points to the time point index; at the start of processing
 *     the line this will be the previous time point index if set previously or
 *     PICKSET_INVALID_INDEX if not; may be modified
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_loadproc_wave(
    PicksetLoadCustomization custom, PicksetLoadCustomWave func, void* arg
);

/*
 * This is a function compatible for use with pickset_add_loadproc_wave.
 * If the string starting with line[*p_offset] has the form
 * 'W %d T %d[whitespace or null]', it will return IVE_ERR_NONE, set
 * *p_wave to the value of the first scanned integer, *p_time to the value
 * of the second scanned integer, and set *p_offset to the index of the
 * first character after the second integer.  Otherwise, it will return
 * IVE_ERR_CORRUPT and leave *p_time, *p_wave, and *p_offset unchanged.
 *
 * arg is not used.
 */
int pickset_simple_start_wave_time(
    int line_number,
    char* line,
    Pickset pickset,
    void* arg,
    int* p_offset,
    int* p_wave,
    int* p_time
);

/*
 * Adds a function to the list of functions that is called by pickset_load
 * when a line beginning with 'L' is encountered.  The functions are called
 * in the order they are registered or until *p_offset is set to something
 * less than zero or to the something greater than or equal to the number of
 * characters in the current line.  pickset_load looks at the integer
 * returned by the last function called.  If it is something other than
 * IVE_ERR_NONE, the load is aborted and that value is returned to the caller
 * of pickset_load.  The function must accept eight arguments:
 *
 *     an int which is the line number currently being processed
 *
 *     a char* which is the NULL terminated string holding the contents of
 *     the current line
 *
 *     a Pickset which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_loadproc_list is called
 *
 *     an int* which points to the offset from the start of the line where
 *     the routine should start processing and should be set before the
 *     function returns to the offset where the next routine should start
 *     processing (or a negative value or value greater than equal to the
 *     length of the line to terminate processing)
 *
 *     an int which is the current wavelength index or PICKSET_INVALID_INDEX
 *     if it has not been set
 *
 *     an int which is the current time point index or PICKSET_INVALID_INDEX
 *     if it has not been set
 *
 *     a PicksetIter* which points to an invalid iterator when the function
 *     is called; if it it still points to an invalid iterator when the
 *     functions returns, the subsequent points (up to the next line
 *     beginning with 'L' or 'W') are read and discarded without invoking the
 *     related callbacks; otherwise, the subsequent points are added to the
 *     list referred to by the iterator
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_loadproc_list(
    PicksetLoadCustomization custom, PicksetLoadCustomList func, void* arg
);

/*
 * This is a function compatible with pickset_add_loadproc_list.  arg
 * is assumed to be NULL or point to a PicksetWTBounds structure.  If
 * arg is NULL no bounds checking is done on iwave and itime.  If arg is
 * non-NULL, the following operations are performed:
 *
 *     if iwave == PICKSET_INVALID_INDEX and arg->flags & PICKSET_WAVE_UNINIT
 *     is non-zero, iwave is set to arg->wave_uninit
 *
 *     if 0 <= iwave + arg->wave_offset < arg->wave_count and
 *     arg->waves_used[iwave + arg->wave_offset] is non-zero, then iwave
 *     is assumed to be valid otherwise it is invalid.
 *
 *     if itime == PICKSET_INVALID_INDEX and arg->flags & PICKSET_TIME_UNINIT
 *     is non-zero, itime is set to arg->time_uninit
 *
 *     if arg->time_min <= itime <= arg->time_max and itime - arg->time_min %
 *     arg->time_step is none-zero, then itime is assumed to be valid otherwise
 *     it is invalid
 *
 *     if iwave or itime is invalid, *p_offset and *p_iter are not modified
 *     and the function returns.  The return valus is IVE_ERR_INVALID if iwave
 *     is invalid and arg->flags & PICKSET_WAVE_OBOUND is zero or if itime is
 *     invalid and arg->flags & PICKSET_TIME_OBOUND is zero; otherwise, it is
 *     IVE_ERR_NONE
 *
 * A list is created for wavelength, iwave, and time point, itime. The created
 * list is assigned the smallest non-negative list number which is not in use;
 * if no non-negative list number is available, *p_offset and *p_iter are not
 * modified and IVE_ERR_OVERFLOW is returned.  If list creation succeeds,
 * *p_iter is set to an iterator to the new list, *p_offset is not modified,
 * and IVE_ERR_NONE is returned; if it fails, *p_offset and *p_iter are not
 * modified and IVE_ERR_NO_MEM is returned.
 */
int pickset_simple_start_list(
    int line_number,
    char* line,
    Pickset pickset,
    void* arg,
    int* p_offset,
    int iwave,
    int itime,
    PicksetIter* p_iter
);

/*
 * Adds a function to the list of functions that is called by pickset_load
 * when loading a line that should specify new point in a list.  The functions
 * are called in the order they are registered until there are no more or
 * until *p_offset is set to something less than zero or to the something
 * greater than or equal to the number of characters in the current line.
 * pickset_load looks at the integer returned by the last function called.
 * If it is something other than IVE_ERR_NONE, the load is aborted and that
 * value is returned to the caller of pickset_load.  The function must accept
 * seven arguments:
 *
 *     an int which is the line number currently being processed
 *
 *     a char* which is the NULL terminated string holding the contents of
 *     the current line
 *
 *     a Pickset which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_loadproc_point is called
 *
 *     an int* which points to the offset from the start of the line where
 *     the routine should start processing and should be set before the
 *     function returns to the offset where the next routine should start
 *     processing (or a negative value or value greater than equal to the
 *     length of the line to terminate processing)
 *
 *     a PicksetIter which refers to the list being read
 *
 *     a PicklistIter* which at the start of processing the line points to
 *     a non-dereferenceable iterator for the point list; may be modified
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_loadproc_point(
    PicksetLoadCustomization custom, PicksetLoadCustomPoint func, void* arg
);

/*
 * This is a function compatible with pickset_add_loadproc_point.  It
 * attempts to scan three floating point values from the string starting at
 * line[*p_offset].  If it is unsucessful, it leaves *p_offset and *p_iter
 * unmodified and returns IVE_ERR_CORRUPT.  If successful, it attempts to
 * create a new point and append it to list.  If this fails, *p_offset is set
 * to -1, *p_iter is unmodified and IVE_ERR_NO_MEM is returned.  Otherwise
 * it sets *p_iter to point to the new point and attempts to scan for an
 * additional integer. If present, the point's name is set to that value,
 * *p_offset is set to the index of the first character after the integer,
 * and IVE_ERR_NONE is returned.  If not present, the point's name is set to
 * zero, *p_offset is set to the index of the first character after the
 * third coordinate, and IVE_ERR_NONE is returned.
 */
int pickset_read_point_coord_name(
    int line_number,
    char* line,
    Pickset pickset,
    void* arg,
    int* p_offset,
    PicksetIter list_iter,
    PicklistIter* p_iter
);

/*
 * Adds a function to the list of functions that is called by pickset_load
 * when loading a line that does not start with '#', 'L', or 'W' and does
 * not appear in a list (i.e. between a line starting with 'W' and a line
 * starting with 'L' or between the start of the file and a line that starts
 * with 'L').  The functions are called in the order they are registered until
 * there are no more or until *p_offset is set to something less than zero or
 * to the something greater than or equal to the number of characters in the
 * current line.  pickset_load looks at the integer returned by the last
 * function called.  If it is something other than IVE_ERR_NONE, the load is
 * aborted and that value is returned to the caller of pickset_load.  The
 * function must accept five arguments:
 *
 *     an int which is the line number currently being processed
 *
 *     a char* which is the NULL terminated string holding the contents of
 *     the current line
 *
 *     a Pickset which is the set of lists being operated upon
 *
 *     a void* which is the value passed as arg when
 *     pickset_set_loadproc_point is called
 *
 *     an int* which points to the offset from the start of the line where
 *     the routine should start processing and should be set before the
 *     function returns to the offset where the next routine should start
 *     processing (or a negative value or value greater than equal to the
 *     length of the line to terminate processing)
 *
 * Returns zero if successful or one if there was insufficient memory to
 * add the function.
 */
int pickset_add_loadproc_misc(
    PicksetLoadCustomization custom, PicksetLoadCustomMisc func, void* arg
);

/*
 * This is a function which is compatible with pickset_add_loadproc_misc.
 * It sets *p_offset to -1 and returns IVE_ERR_CORRUPT.  arg is not used.
 */
int pickset_misc_is_error(
    int line_number, char* line, Pickset pickset, void* arg, int* p_offset
);

/*
 * Allocates space for a property to be associated with a Pickset.
 * If there was other data with the same property code as property_code, that
 * data is overwritten.  type_code may be used by applications to indicate the
 * data type used for the property.
 *
 * Returns a pointer to the unitialized property data if successful or NULL
 * if not.
 */
void* pickset_allocate_property(
    Pickset pickset,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size
);

/*
 * Copies n_elements * size_t bytes from data into property, property_code,
 * associated with the Picklist, list.  Any exising data for that property is
 * overwritten.  Returns zero if successful or 1 if memory for the property
 * could not be allocated.
 */
int pickset_set_property(
    Pickset pickset,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size,
    const void* data
);

/*
 * Removes the given property, if it exists, from the Pickset, pickset.
 */
void pickset_delete_property(Pickset pickset, int property_code);

/*
 * Returns a pointer to property data matching the given property_code for
 * the Pickset, pickset.  A non-NULL value is returned if the property exists,
 * and the values pointed to by p_type_code, p_element_size, p_n_elements,
 * if non-NULL, will be set to the property's attributes.  The returned
 * pointer may be used to modify the property.  If the list did not have
 * the property defined, NULL is returned and the values pointed to by
 * p_type_code, p_element_size, and p_n_elements are not modified.
 */
void* pickset_get_property(
    Pickset pickset,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

const void* pickset_get_property_ro(
    PicksetRO pickset,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

/*
 * Attempts to add a new list with for the given wavelength and time point
 * and with the given list number.  Returns a valid iterator if successful.
 * Otherwise an invalid iterator is returned and if p_reason is non_NULL
 * sets *p_reason to one of the following:
 *
 * IVE_ERR_NO_MEM
 *      There was insufficient space or system resources.
 * IVE_ERR_EXIST
 *      A list with the same wavelength, time, and list indices already exists.
 * IVE_ERR_INVALID
 *      wave_index, time_index, or list_index is equal to
 *      PICKSET_INVALID_INDEX.
 *
 * The behavior is undefined if pickset is not valid.
 */
PicksetIter pickset_insert(
    Pickset pickset,
    int wave_index,
    int time_index,
    int list_index,
    int* p_reason
);

/*
 * Erases the list, if any, with the wavelength, time, and list indices.
 */
void pickset_erase(
    Pickset pickset, int wave_index, int time_index, int list_index
);

/*
 * Returns an iterator to the list with the given wavelength, time, and
 * list indices or returns an invalid iterator if no such list is available.
 * The results are undefined if pickset is not valid.
 */
PicksetIter pickset_search(
    Pickset pickset, int wave_index, int time_index, int list_index
);
PicksetIterRO pickset_search_ro(
    PicksetRO pickset, int wave_index, int time_index, int list_index
);

/*
 * Returns an iterator to the first (in the default ordering by the wavelength,
 * time, and list indices) list in the set or an invalid iterator if the
 * set is empty.  The results are undefined if pickset is not valid.
 */
PicksetIter pickset_first(Pickset pickset);
PicksetIterRO pickset_first_ro(PicksetRO pickset);

/*
 * Returns an iterator to the last (in the default ordering by the wavelength,
 * time, and list indices) list in the set or an invalid iterator if the
 * set is empty.  The results are undefined if pickset is not valid.
 */
PicksetIter pickset_last(Pickset pickset);
PicksetIterRO pickset_last_ro(PicksetRO pickset);

/*
 * Returns an invalid iterator for the given Pickset.  The results are
 * undefined if pickset is not valid.
 */
PicksetIter pickset_invalid_iter(Pickset pickset);
PicksetIterRO pickset_invalid_iter_ro(PicksetRO pickset);

/*
 * Returns true if the given iterator is dereferenceable or false if not.
 */
int pickset_iter_valid(PicksetIter iter);
int pickset_iter_valid_ro(PicksetIterRO iter);

/*
 * Returns true if the two iterators point to the same list and false if
 * they do not.
 */
int pickset_iter_equal(PicksetIter iter1, PicksetIter iter2);
int pickset_iter_equal_vro(PicksetIter iter1, PicksetIterRO iter2);
int pickset_iter_equal_rov(PicksetIterRO iter1, PicksetIter iter2);
int pickset_iter_equal_roro(PicksetIterRO iter1, PicksetIterRO iter2);

/*
 * Converts a read/write iterator to a read-only one.
 */
PicksetIterRO pickset_iter_to_ro(PicksetIter iter);

/*
 * Returns the list pointed to by the given iterator.  The behavior is
 * undefined if the iterator is not dereferenceable.
 */
Picklist pickset_iter_dereference(PicksetIter iter);
PicklistRO pickset_iter_dereference_ro(PicksetIterRO iter);

/*
 * Returns an iterator to the next list in the default ordering.  If the
 * iterator is invalid, returns the equivalent of pickset_first().
 */
PicksetIter pickset_iter_increment(PicksetIter iter);
PicksetIterRO pickset_iter_increment_ro(PicksetIterRO iter);

/*
 * Returns an iterator to the previous list in the default ordering.  If the
 * iterator is invalid, returns the equivalent of pickset_last().
 */
PicksetIter pickset_iter_decrement(PicksetIter iter);
PicksetIterRO pickset_iter_decrement_ro(PicksetIterRO iter);

/*
 * Removes the list pointed to by the iterator, iter.  The behavior is
 * undefined if the iterator is not dereferenceable.
 */
void pickset_iter_remove(PicksetIter iter);

/*
 * Returns the number of points in a list.
 */
int picklist_size(PicklistRO list);

/*
 * Sets *p_wave_index to the wavelength index for the list, *p_time_index
 * to the time index for the list, and *p_list_index to the list number for
 * the list.
 */
void picklist_indices(
    PicklistRO list, int* p_wave_index, int* p_time_index, int* p_list_index
);

/*
 * Removes all points in a Picklist.  Invalidates any iterators pointing
 * to the list.
 */
void picklist_clear(Picklist list);

/*
 * Allocates space for a property to be associated with a Picklist.
 * If there was other data with the same property code as property_code, that
 * data is overwritten.  type_code may be used by applications to indicate the
 * data type used for the property.
 *
 * Returns a pointer to the unitialized property data if successful or NULL
 * if not.
 */
void* picklist_allocate_property(
    Picklist list,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size
);

/*
 * Copies n_elements * size_t bytes from data into property, property_code,
 * associated with the Picklist, list.  Any exising data for that property is
 * overwritten.  Returns zero if successful or 1 if memory for the property
 * could not be allocated.
 */
int picklist_set_property(
    Picklist list,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size,
    const void* data
);

/*
 * Removes the given property, if it exists, from the Picklist, list.
 */
void picklist_delete_property(Picklist list, int property_code);

/*
 * Returns a pointer to property data matching the given property_code for
 * the Picklist, list.  A non-NULL value is returned if the property exists,
 * and the values pointed to by p_type_code, p_element_size, p_n_elements,
 * if non-NULL, will be set to the property's attributes.  The returned
 * pointer may be used to modify the property.  If the list did not have
 * the property defined, NULL is returned and the values pointed to by
 * p_type_code, p_element_size, and p_n_elements are not modified.
 */
void* picklist_get_property(
    PicklistRO list,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

const void* picklist_get_property_ro(
    Picklist list,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

/*
 * Returns an iterator to the points in a Picklist.
 */
PicklistIter picklist_start(Picklist list);
PicklistIterRO picklist_start_ro(PicklistRO list);

/*
 * Returns an iterator one past the end of the points in a Picklist.
 */
PicklistIter picklist_end(Picklist list);
PicklistIterRO picklist_end_ro(PicklistRO list);

/*
 * Adds count points (whose coordinates and names are unitialized and
 * initially have no additional properties), at the end of the given
 * list.  Returns an iterator to the first of the inserted points if
 * successful or a non-dereferenceable iterator if not.  If successful
 * invalidates any other iterators to the list.
 */
PicklistIter picklist_append(Picklist list, int count);

/*
 * Adds count points (whose coordinates and names are unitialized and
 * initially have no additional properties), at the beginning of the given
 * list.  Returns an iterator to the first of the inserted points if
 * successful or a non-dereferenceable iterator if not.  If successful
 * invalidates any other iterators to the list.
 */
PicklistIter picklist_prepend(Picklist list, int count);

/*
 * Returns true if the given iterator is dereferenceable and false if not.
 */
int picklist_iter_valid(PicklistIter iter);
int picklist_iter_valid_ro(PicklistIterRO iter);

/*
 * Returns true if the two iterators point to the same position and false if
 * they do not.
 */
int picklist_iter_equal(PicklistIter iter1, PicklistIter iter2);
int picklist_iter_equal_vro(PicklistIter iter1, PicklistIterRO iter2);
int picklist_iter_equal_rov(PicklistIterRO iter1, PicklistIter iter2);
int picklist_iter_equal_roro(PicklistIterRO iter1, PicklistIterRO iter2);

/*
 * Converts a read/write iterator to a read-only one.
 */
PicklistIterRO picklist_iter_to_ro(PicklistIter iter);

/*
 * Returns the a pointer to the point pointed to by iter.  The behavior is
 * undefined if the iterator is not dereferenceable.
 */
Pickpt* picklist_iter_dereference(PicklistIter iter);
const Pickpt* picklist_iter_dereference_ro(PicklistIterRO iter);

/*
 * Inserts count points (whose coordinates and names are unitialized and
 * initially have no additional properties), before the point pointed to
 * by iter.  Returns an iterator to the first of the inserted points if
 * successful or a non-dereferenceable iterator if not.  If successful,
 * invalidates any iterators to the list.
 */
PicklistIter picklist_iter_insert(PicklistIter iter, int count);

/*
 * Swaps the point referenced by iter_a with the point referenced by iter_b.
 * The result is undefined if either of the iterators is not dereferenceable
 * or do not refer to the same list.
 */
void picklist_iter_swap(PicklistIter iter_a, PicklistIter iter_b);

/*
 * Deletes the point pointed to by iter.  Invalidates any iterators
 * to the list.
 */
void picklist_iter_erase(PicklistIter iter);

/*
 * Deletes all points starting with the one pointed to by start up to but
 * not including the one pointed to by end.  Invalidates any iterators
 * to the list.
 */
void picklist_iter_erase_range(PicklistIter start, PicklistIter end);

/*
 * Returns an iterator pointing to the next point.
 */
PicklistIter picklist_iter_increment(PicklistIter iter);
PicklistIterRO picklist_iter_increment_ro(PicklistIterRO iter);

/*
 * Returns an iterator pointing to the previous point.
 */
PicklistIter picklist_iter_decrement(PicklistIter iter);
PicklistIterRO picklist_iter_decrement_ro(PicklistIterRO iter);

/*
 * Records the fact that the point is drawn in the image window whose number
 * is window_number with the graphic whose ID is id.  Returns 0 if successful
 * and 1 if not.  The image window display is not updated.
 */
int picklist_iter_set_graphic_id(PicklistIter iter, int id, int window_number);

/*
 * Removes the graphic, if any, associated with the point in the image window
 * whose number is window_number.  If window_number is IW_ALL_WINDOWS, the
 * graphics for the point in all windows are removed.  The image window
 * display is not updated.
 */
void picklist_iter_remove_graphic(PicklistIter iter, int window_number);

/*
 * Allocates space for a property to be associated with the point pointed to
 * by iter.  If there was other data with the same property code as
 * property_code, that data is overwritten.  type_code may be used by
 * applications to indicate the data type used for the property.
 *
 * Returns a pointer to the unitialized property data if successful or NULL
 * if not.
 */
void* picklist_iter_allocate_property(
    PicklistIter iter,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size
);

/*
 * Copies n_elements * size_t bytes from data into property, property_code,
 * associated with the point pointed to by iter.  Any exising data for that
 * property is overwritten.  Returns zero if successful or 1 if memory for
 * the property could not be allocated.
 */
int picklist_iter_set_property(
    PicklistIter list,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size,
    const void* data
);

/*
 * Removes the given property, if it exists, from the point pointed to by iter.
 */
void picklist_iter_delete_property(PicklistIter iter, int property_code);

/*
 * Returns a pointer to property data matching the given property_code for
 * the point pointed to by iter.  A non-NULL value is returned if the property
 * exists, and the values pointed to by p_type_code, p_element_size,
 * p_n_elements, if non-NULL, will be set to the property's attributes.  The
 * returned pointer may be used to modify the property.  If the point did not
 * have the property defined, NULL is returned and the values pointed to by
 * p_type_code, p_element_size, and p_n_elements are not modified.
 */
void* picklist_iter_get_property(
    PicklistIter list,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

const void* picklist_iter_get_property_ro(
    PicklistIterRO iter,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

#ifdef __cplusplus
}
#endif

#endif
