#ifndef UCSF_MSG_IWOBJMONINTERFACE_H
#define UCSF_MSG_IWOBJMONINTERFACE_H

/*
 * Declares a set of functions that perform the basic manipulations
 * that Priism's Monitor needs so that the Monitor does not directly
 * use the C++ internals of the object library.
 */

#include "ive_shm.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 * header_offsets is the 0-terminated array of shared memory offsets for
 * the IWGrHeader data structures.  arena is the shared memory arena
 * holding those structures.  wid is the Monitor's window number.
 * zwt holds the z, wavlength, and time indices for the section to be
 * drawn.  mode should be IWOBJ_REAL_SP or IWOBJ_DATA_SP to select which
 * objects are to be drawn.  color_func is a function pointer to be invoked to
 * change the rendering colors for the image window.
 */
void IWMonDrawObjectsForSection(
    IVEArena arena,
    IVEPtrRep offset_to_header_list,
    int wid,
    const int zwt[3],
    int mode,
    void (*color_func)(int wid, int color_index)
);

/*
 * p_offset_to_header_list is a pointer to the IVEPtrRep holding the shared
 * memory offset for the start of the list of offsets to the header data
 * structures.  arena is the shared memory arena relative to which all the
 * offsets are calculated.  wid is the Monitor's window number.
 */
void IWMonDetachAllObjects(
    IVEArena arena, IVEPtrRep* p_offset_to_header_list, int wid
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
