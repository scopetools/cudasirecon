#ifndef UCSF_MSG_POLYFUNCS_H
#define UCSF_MSG_POLYFUNCS_H

#include "IWInclude.h"  /* IW_POINT */


/*
 * A Polylist is a list of polygons by section number.  It is a associated
 * with an IM library stream.  The Polylist may be private (only visible
 * to the application that created it) or shared.  Shared lists must
 * be associated with a stream attached to an image window.  After the
 * Polylist is created, the Polylist calls do not monitor if the attributes
 * of the associated stream change or if it is closed:  this is the
 * responsibility of the application.
 *
 * To access the polygons in the list, the iterators, Polyiter or
 * PolyiterConst, are used.  Functions are defined for these iterators
 * to modify or retrieve the polygon attributes,  insert or delete elements
 * into the polygon list, and handle the graphics for the polygons.
 * Additionally there are functions so applications may associate
 * application-specific data, i.e. a property, with a polygon.  Properties
 * are referred to by an integer code and an integer code is also associated
 * with the data type used for the property.  Certain predefined properties
 * types and the corresponding data structures are defined here (look for the
 * POLYPROP_* and POLYTYPE_* constants).
 *
 * Some utility functions to work with the vertices for a polygon are
 * provided: poly_compute_bounding_box, poly_is_inside,
 * poly_determine_inner_pts.
 *
 */


struct PolylistImpl;
typedef struct PolylistImpl* Polylist;
 
struct PolyImpl;
typedef struct Polyiter {
    struct PolyImpl* p_poly;
    Polylist list;
    int section;
} Polyiter;
typedef struct PolyiterConst {
    const struct PolyImpl* p_poly;
    const struct PolylistImpl* list;
    int section;
} PolyiterConst;

/*
 * Specifies the bounding box for the polygon (in the same coordinates
 * which the vertices are specified.  The associated data should be one
 * PolyBoundingBox (with a type code POLYTYPE_BOUNDING_BOX).
 */
#define POLYPROP_BOUNDING_BOX 1
/*
 * Specifies the starting point that was used to automatically generate
 * the polygon.  The associated data should be on PolyPt (with a type code
 * of POLYTYPE_PT).
 */
#define POLYPROP_SEED_PT 2
/*
 * Specifies the one dimensional coordinates for the pixels in the
 * interior of the polygon.  The associated data will be zero or more
 * integers (with a type code POLYTYPE_COORD_1D).
 */
#define POLYPROP_INTERNAL_PTS 3

#define POLYTYPE_BOUNDING_BOX 1
#define POLYTYPE_PT 2
#define POLYTYPE_COORD_1D 3

typedef struct PolyBoundingBox {
    float x_min;
    float x_max;
    float y_min;
    float y_max;
} PolyBoundingBox;

typedef struct PolyPt {
    int x;
    int y;
} PolyPt;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Creates a new polygon list or attaches to an existing one.  istream must
 * be an open IM stream for which IMRdHdr has been called.  flags must
 * be a bitwise or of one or more of the following: 0, POLYLIST_PRIVATE,
 * POLYLIST_SHARED, POLYLIST_CREATE, POLYLIST_EXCL, POLYLIST_CLEAR,
 * POLYLIST_WLOCK.  Normally, whether a list is private or shared is
 * determined from the stream number.  To force a list to be treated as
 * private set POLYLIST_PRIVATE in flags and to force it to be shared set
 * POLYLIST_SHARED.  If POLYLIST_CREATE is set, a new list will be created if
 * necessary.  Since no provision is made for attaching to an existing
 * private list, POLYLIST_CREATE should be specified when creating private
 * lists.  If POLYLIST_EXCL is set and a shared list is to be created,
 * the function call will fail if a list is already associated with the
 * window.  If POLYLIST_CLEAR is set and the call attaches to an existing
 * shared list, the list will be cleared before polylist_attach returns.
 * If POLYLIST_WLOCK is set, the list will be locked for writing when
 * polylist_attach returns; it is the caller's responsibility to release the
 * lock when it is no longer needed.
 *
 * Returns a non-NULL value if the function call succeeded and a NULL
 * value if it failed.  If the application wants to determine the cause
 * of failure, it should pass a non-NULL value for p_status.  Then, if
 * the call fails, *p_status will be set to one of the following:
 *
 * IVE_ERR_NO_MEM
 *      There was insufficient space or system resources.
 * IVE_ERR_EXIST
 *      POLYLIST_EXCL was set but a list with the same identifying
 *      characteristics was found.
 * IVE_ERR_INVALID
 *      POLYLIST_CREATE was not set and there was no list to attach to,
 *      POLYLIST_SHARED was specified and istream was not attached to a window,
 *      or POLYLIST_SHARED and POLYLIST_PRIVATE were specified.
 *
 * polylist_attach may deadlock if the caller hold a read or write lock
 * to a shared list which is asoociated with the same window as istream.
 */
enum {
    POLYLIST_PRIVATE = 1,
    POLYLIST_SHARED = 2,
    POLYLIST_CREATE = 4,
    POLYLIST_EXCL = 8,
    POLYLIST_CLEAR = 16,
    POLYLIST_WLOCK = 32
};
Polylist polylist_attach(int flags, int istream, int* p_status);

/*
 * Detaches from a polygon list.  If this is the last application attached
 * to the list, the resources associated with the list are released.
 *
 * polylist_detach may not be called if the caller has locked the list.
 */
void polylist_detach(Polylist list);

/*
 * Deletes the resources for a polygon list.  If the list is shared, this
 * may be dangerous (there may be other applications still using the list).
 *
 * polylist_delete may not be called if the caller has locked the list.
 */
void polylist_delete(Polylist list);

/*
 * Lock or unlock the list for reading or writing.  These calls have no
 * effect on  private list.
 */
void polylist_lock_read(Polylist list);
void polylist_lock_write(Polylist list);
void polylist_unlock_read(Polylist list);
void polylist_unlock_write(Polylist list);

/*
 * Returns true if the list is shared, false if it is not.
 */
int polylist_shared(const Polylist list);

/*
 * Returns the stream number associated with the list or 0 if list is invalid.
 */
int polylist_stream(const Polylist list);

/*
 * Saves the polygon list to a file.  Returns IVE_ERR_NONE if successful
 * or one of the following if not:
 *
 *      IVE_ERR_BAD_FILENAME
 *           file_name could not be opened for writing
 *      IVE_ERR_IO
 *           An I/O error occurred while writing the file.
 *      IVE_ERR_INVALID
 *           list was invalid.
 */
int polylist_save(const Polylist list, const char* file_name);

/*
 * Loads a polygon list from a file.  If append is true, the list will
 * be added to what is already present, otherwise the current contents
 * are cleared before the list is loaded.  Returns IVE_ERR_NONE if successful
 * one of the following if not:
 *
 *      IVE_ERR_BAD_FILENAME
 *           file_name could not be opened for reading.
 *      IVE_ERR_IO
 *           An I/O error occurred while reading the file.
 *      IVE_ERR_CORRUPT
 *           The polygon file did not have the correct format.
 *      IVE_ERR_NO_MEM
 *           There was insufficient space to store the loaded polygons.
 *      IVE_ERR_INVALID
 *           list was invalid.
 *
 * If the function call is not successful the contents of the list are not
 * modified.
 */
int polylist_load(Polylist list, const char* file_name, int append);

/*
 * The polygons in the list are stored by section number.  This
 * section number may not be the same as the corresponding section number
 * for the IM stream.  polylist_zwt_to_index converts z, wavelength,
 * and time indices to a corresponding section number usable with the
 * polylist calls.  -1 is returned if the list is invalid or the given
 * indices are out of bounds for the list.  polylist_index_to_zwt returns the
 * z, wavelength, and time indices corresponding to a section number used by
 * the polylist calls.  The indices will be set to -1 if the list is invalid
 * or the section index is out of bounds.
 */
int polylist_zwt_to_index(const Polylist list, int z, int w, int t);
void polylist_index_to_zwt(
    const Polylist list, int index, int* p_z, int* p_w, int* p_t
);

/*
 * Returns the total number of sections in the polygon list.
 */
int polylist_section_count(const Polylist list);

/*
 * Return an iterator to the polygons in a single section.
 */
Polyiter polylist_iterator(Polylist list, int section);
PolyiterConst polylist_ro_iterator(const Polylist list, int section);

/*
 * Return an iterator point to one past the end of the polygons in a single
 * section.
 */
Polyiter polylist_iterator_end(Polylist list, int section);
PolyiterConst polylist_ro_iterator_end(const Polylist list, int section);

/*
 * Returns true if the iterator, iter, is past the end of the list or
 * false if it is not.
 */
int polyiter_at_end(Polyiter iter);
int polyiter_ro_at_end(PolyiterConst iter);

/*
 * Converts a read/write iterator to a read-only one.
 */
PolyiterConst polyiter_to_ro(Polyiter iter);

/*
 * Returns true if the two iterators point to the same position in
 * a polygon list and false if they do not.
 */
int polyiter_equal(Polyiter iter1, Polyiter iter2);
int polyiter_vro_equal(Polyiter iter1, PolyiterConst iter2);
int polyiter_rov_equal(PolyiterConst iter1, Polyiter iter2);
int polyiter_roro_equal(PolyiterConst iter1, PolyiterConst iter2);

/*
 * Retrieve or modify the attributes of the polygon pointed to by
 * iter.  For all except, *_section_index, the results are undefined
 * if the functions are called with an iterator for which
 * polyiter_at_end is true.
 */
int polyiter_section_index(Polyiter iter);
int polyiter_ro_section_index(PolyiterConst iter);

int polyiter_closed(Polyiter iter);
int polyiter_ro_closed(PolyiterConst iter);
void polyiter_set_closed(Polyiter iter, int closed);

int polyiter_allwaves(Polyiter iter);
int polyiter_ro_allwaves(PolyiterConst iter);
void polyiter_set_allwaves(Polyiter iter, int all);

int polyiter_vertex_count(Polyiter iter);
int polyiter_ro_vertex_count(PolyiterConst iter);

IW_POINT* polyiter_vertices(Polyiter iter);
const IW_POINT* polyiter_ro_vertices(PolyiterConst iter);

/*
 * Changes the number of vertices to count and returns a pointer
 * to the start of the vertex array if successful or a NULL if it was
 * not successful or the new count was zero.  The contents of the vertex
 * array up to the minimum of the previous vertex count and the new vertex
 * count are preserved after a successful call, the vertices are unmodified
 * after an unsuccessful one.  Expanding the vertex array may invalidate any
 * pointer previously returned by polyiter_vertices, polyiter_ro_vertices,
 * polyiter_set_vertex_count, or polyiter_insert_vertices for this polygon.
 * The coordinates for any added vertices are not initialized. The result
 * is undefined if polyiter_at_end would return true for iter.
 */
IW_POINT* polyiter_set_vertex_count(Polyiter iter, int count);

/*
 * Deletes count vertices starting at the given position from the polygon
 * referred to by iter.  The result is undefined if polyiter_at_end
 * would return true for iter.
 */
void polyiter_delete_vertices(Polyiter iter, int position, int count);

/*
 * Creates space for count vertices starting at the given position.  The
 * original vertices at those positions are shifted.  Returns a pointer
 * to the start of the vertex array if successful or a NULL if not.
 * The contents of the vertex array array are not changed if the call is
 * not successful.  May invalidate any pointer previously returned by
 * polyiter_vertices, polyiter_ro_vertices, polyiter_set_vertex_count,
 * or polyiter_insert_vertices for this polygon. The coordinates for any
 * added vertices are not initialized. The result is undefined if
 * polyiter_at_end would return true for iter.
 */
IW_POINT* polyiter_insert_vertices(Polyiter iter, int position, int count);

/*
 * Allocates space for a property to be associated with the polygon pointed to
 * by iter.  If there was other data with the same property code as
 * property_code, that data is overwritten.  type_code may be used by
 * applications to indicate the data type used for the property.
 *
 * Returns a pointer to the unitialized property data if successful or NULL
 * if not.  The result is undefined if polyiter_at_end would return true for
 * iter.
 */
void* polyiter_allocate_property(
    Polyiter iter,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size
);

/*
 * Copies n_elements * size_t bytes from data into property, property_code,
 * associated with the polygon pointed to by iter.  Any exising data for
 * that property is overwritten.  Returns zero if successful or 1 if
 * memory for the property could not be allocated.  The result is undefined
 * if polyiter_at_end would return true for iter.
 */
int polyiter_set_property(
    Polyiter iter,
    int property_code,
    int type_code,
    size_t n_elements,
    size_t element_size,
    const void* data
);

/*
 * Removes the given property, if it exists, from the polygon pointed to
 * by iter.  The result is undefined if polyiter_at_end would return true
 * for iter.
 */
void polyiter_delete_property(Polyiter iter, int property_code);

/*
 * Returns a pointer to property data matching the given property_code for
 * the polygon pointed to by iter.  A non-NULL value is returned if the
 * property exists, and the values pointed to by p_type_code, p_element_size,
 * p_n_elements, if non-NULL, will be set to the property's attributes.
 * The returned pointer may be used to modify the property.  If the polygon
 * did not have the property defined, NULL is returned and the values
 * pointed to by p_type_code, p_element_size, and p_n_elements are not
 * modified.
 *
 * The result is undefined if polyiter_at_end would return true for iter.
 */
void* polyiter_get_property(
    Polyiter iter,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

const void* polyiter_ro_get_property(
    PolyiterConst iter,
    int property_code,
    int* p_type_code,
    size_t* p_n_elements,
    size_t* p_element_size
);

/*
 * Advances the iterator to the next position and returns the a copy
 * of the incremented interator.  The result is undefined
 * if polyiter_at_end would return true for iter.
 */
Polyiter polyiter_increment(Polyiter iter);
PolyiterConst polyiter_ro_increment(PolyiterConst iter);

/*
 * Deletes the polygon currently pointed to by iter and advances the
 * iter to point to the next polygon.  Any graphics associated with
 * the polygon are removed but the window is not redrawn.  A copy of
 * the incremented iterator is returned.  After the call, operations on
 * iterators which point to the deleted polygon will give undefined results.
 */
Polyiter polyiter_remove_increment(Polyiter iter);

/*
 * Attempts to create a new polygon.  If this succeeds, the newly
 * created polygon is inserted immediately before the position pointed
 * to by iter, iter is modified to point to the new polygon, and
 * zero is returned.  Otherwise, iter is not changed and a non-zero
 * value is returned.  A newly created polygon has zero vertices, is
 * closed, and does not appear in all waves.
 */
int polyiter_insert(Polyiter* p_iter);

/*
 * Overwrites the polygon pointed to by iter2 by the polygon pointed
 * to by iter1.  If iter2 pointed to a polygon with graphics, those
 * graphics are removed but no graphics corresponding to the new shape
 * are added and the window is not redrawn.  The properties, if any,
 * associated with the polygon pointed to by iter1 are not copied.
 * The properties associated with the polygon pointed to by iter2 are
 * removed.  Returns zero if successful and a non-zero value if not.  If
 * the call fails, the polygon pointed to by iter2 is not modified.
 */
int polyiter_copy(Polyiter iter1, Polyiter iter2);
int polyiter_ro_copy(PolyiterConst iter1, Polyiter iter);

/*
 * For the polygon currently pointed to by iter, removes any associated
 * graphics and adds new graphics corresponding to its current shape.
 * If all_zwt is true, the graphics will be added to all sections.
 * If selected_vertex is a valid vertex number (between 0 and
 * polyiter_vertex_count(iter) - 1) then that vertex will be drawn with an X.
 * The image window is not redrawn.  The result is undefined if
 * polyiter_at_end would return true for iter.
 */
void polyiter_draw(Polyiter iter, int color, int all_zwt, int selected_vertex);

/*
 * If the polygon pointed to by iter was previously drawn, cause the
 * graphics to be displayed on the image window.
 */
void polyiter_display(Polyiter iter);

/*
 * For the polygon currently pointed to by iter, removes any associated
 * graphics.  The image window is not redrawn.  The result is undefined
 * if polyiter_at_end would return true for iter.
 */
void polyiter_hide(Polyiter iter);

/*
 * Preconditions
 *      pts contains vertex_count points which all lie in the same z plane.
 *
 * Postconditions
 *      Returns false if vertex_count is less than 1.  If vertex_count is
 *      one, returns true if include_edge is true and (test_pt[0], test_pt[1])
 *      is the same as (pts[0].x, pts[0].y); otherwise returns false.  If
 *      vertex_count is two, returns true if include_edge is true and
 *      (test_pt[0], test_pt[1]) falls on the line segment between
 *      (pts[0].x, pts[0].y) and (pts[1].x, pts[1].y).  If vertex_count is
 *      greater than 2, let p be the polygon formed by taking the vertex_count
 *      points from pts in order as the vertices.  Then true is returned if
 *      include_edge is true and the point falls on the boundary of the polygon
 *      or if the point is inside the polygon; otherwise, false is returned.
 */
int poly_is_inside(
    const IW_POINT* pts,
    int vertex_count,
    const float test_pt[3],
    int include_edge
);

/*
 * Preconditions
 *      pts contains vertex_count points.
 *
 * Postconditions:
 *      If vertex_count is less than 1, *p_minx, *p_miny, *p_maxx, *p_maxy
 *      are unchanged.  If vertex_count is greater than or equal to one,
 *      *p_minx is set to the largest possible value for which pts[i].x >=
 *      *p_minx is true for all i in [0, vertex_count - 1].  *p_miny is
 *      set to the largest possible value for which pts[i].y >=
 *      *p_miny is true for all i in [0, vertex_count - 1].  *p_maxx is
 *      set to the smallest possible value for which pts[i].x <= *p_maxx
 *      is true for all i in [0, vertex_count - 1].  *p_maxy is
 *      set to the smallest possible value for which pts[i].y <= *p_maxy
 *      is true for all i in [0, vertex_count - 1].
 */
void
poly_compute_bounding_box(
    const IW_POINT* pts,
    int vertex_count,
    float* p_minx,
    float* p_miny,
    float* p_maxx,
    float* p_maxy
);

/*
 * Determines the one-dimensional coordinates for the pixels that fall
 * within the polygon whose vertices are passed in pts.  If the return
 * value is non-NULL, *p_n_inner is set to the number of interior points
 * and the returned points points to the *p_n_inner one-dimensional
 * coordinates (computed as the x pixel index + y pixel index * width).
 * If allocator is NULL, the memory for the interior points is allocated
 * with malloc() and should be freed by the caller with free(); otherwise
 * the memory is allocated by a call to allocator.  The return value is NULL
 * if the memory allocation failed.
 */
int* poly_determine_inner_points(
    const IW_POINT* pts,
    int vertex_count,
    int width,
    int* p_n_inner,
    void* (*allocator)(size_t n_elements, size_t element_size, void* p_param),
    void* allocator_param
);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
