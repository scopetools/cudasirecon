#ifndef UCSF_MSG_SELECTREGION_H
#define UCSF_MSG_SELECTREGION_H

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Initiates interactive selection of an XY region from an image window.
 *
 * Preconditions:
 *     istr is the number of a stream attached to an image window.  The
 *         header for that stream must have been read in and IW events
 *         enabled on the stream.
 *     box is a four element vector which is interpreted as follows:
 *         (box[0], box[1]) is the lower lefthand corner of the initial
 *         selection region and (box[2], box[3]) is the upper righhand
 *         corner of the selection region.  If box is NULL or the selection
 *         region is invalid (box[0] > box[2] or box[1] > box[3]), the initial
 *         box will be the complete x/y range for the stream.  The limits of
 *         the initial box will be coerced to lie within the x/y range for
 *         the stream.
 *     update_cb is NULL or a function to be invoked whenever the selected
 *         region changes during the selection process.
 *     complete_cb is NULL or a function to be invoked when the mouse
 *         pointer leaves the image window.
 *
 * Postcondition
 *     The function returns a non-zero value if there were insufficient
 *     resources.
 *
 *     Otherwise zero is returned; undefined behavior will result if
 *     istr is closed prior to the region selection process terminating.
 *     The selection process terminates at the earliest of the following
 *     events:
 *          1) the region selection completes and complete_cb is called with
 *             istr as the stream number
 *          2) CancelSelectRegion with istr as the stream number and returns
 *     Calling SelectRegion again before one of the above events occurs
 *     effectively cancels the current region selection process and restarts
 *     it with the given arguments.
 *
 *     After SelectRegion returns zero, update_cb is called while the mouse
 *     pointer is in the image window attached to istr, and the region
 *     selection process has not terminated.  update_cb is passed the stream
 *     number passed to SelectRegion, the bounds of the selected region
 *     ((box[0], box[1]) is the lower lefthand corner and (box[2], box[3]) is
 *     the upper righthand corner), and the value of update_closure passed to
 *     SelectRegion.
 *
 *     complete_cb is called after SelectRegion returns zero and the region
 *     selection process terminates.  The region selection process terminates
 *     at the first occurence of one of the following:
 *
 *          CancelSelectRegion is called on istr
 *
 *          SelectRegion is called again on istr
 *
 *          The mouse button has been pressed while the pointer is in the
 *          image window and subsequently the mouse button was release and
 *          the pointer left the window.
 *
 *     complete_cb is passed the stream number passed to SelectRegion, the
 *     bounds of the selected region ((box[0], box[1]) is the lower lefthand
 *     corner and (box[2], box[3]) is the upper righthand corner, the ID of
 *     the graphic used to mark the selected region boundary, and the value
 *     of complete_clousre passed to SelectRegion.  Once complete_cb is
 *     called, no further operations are called internally on the graphic ID:
 *     i.e. if it is to be removed the caller of SelectRegion must arrange
 *     for that to happen.
 *
 * Side effects:
 *     Because WMCancelEventHandler is used internally, undefined behavior
 *     results if the application attempts to monitor EnterWindow, LeaveWindow,
 *     ButtonPress, ButtonRelease, or PointerMotion events on istr between
 *     the time of the call to SelectRegion and the termination of the
 *     region selection.
 */
int SelectRegion(
    int istr,
    const int box[4],
    void (*update_cb)(int istr, int box[4], void* closure),
    void* update_closure,
    void (*complete_cb)(int istr, int box[4], int graphics_id, void* closure),
    void* complete_closure
);

/*
 * Cancels a region selection initiated by SelectRegion.
 *
 * Postcondition:
 *     The resources associated with region selection are freed and
 *     the callbacks registered with SelectRegion for istr will not be invoked
 *     until a subsequent SelectRegion call reregisters them.  If
 *     window_closed is true, it is assumed that the image window was closed
 *     and that the resources managed by the image window do not need to be
 *     released.
 */
void CancelSelectRegion(int istr, int window_closed);

/*
 * Changes the wavelength used for coordinate calculations during region
 * selection.
 *
 * Precondition:
 *     SelectRegion has been called for istr and returned zero and
 *     the region selection process has not terminated.
 *
 *     wave is an element of [0, number of wavelengths - 1] where the number
 *     of wavelengths in the stream istr.
 *
 * Postcondition:
 *     Coordinate conversions for region selection (to determine the position
 *     of the mouse in the image) are done with respect to the wavelength with
 *     index wave.   By default, coordinate conversions are done with respect
 *     to the lowest index wavelength mapped at the time SelectRegion was
 *     called for istr or the zero index wavelength if no wavelengths were
 *     mapped.
 */
void ChangeSelectRegionWave(int istr, int wave);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
