/*
 * Displays a trace as the mouse is moved across an image window.  The graphics
 * are done with the overlay planes, if available, or the IWBgnFastGr set
 * of calls.  Does not attempt to handle the case where the user changes the
 * current section while the mouse is over the window.
 */

#include "WMInclude.h"
#include "IWInclude.h"
#include <stdio.h>      /* sprintf */
#include <stdlib.h>     /* exit */


static int quit(int* p_unused);
static int EventHandler(int itype, XEvent* event);


static IW_POINT_2D gfPt[1000];
static int giNPt = 0;
static int giGID = -1;


int main(int argc, char* argv[])
{
    static char message[30] = "Displays a trace as the mouse";
    static char message1[30] = "is moved over an image window";

    WMInit("MotionEvent Demo");
    WMSetLoc(500, 500);
    WMAddText(message, sizeof(message) - 1, 1);
    WMNewRow();
    WMAddText(message1, sizeof(message1) - 1, 1); 
    WMNewRow();
    WMAddFuncButton("Exit", quit, NULL, 0, 0);
    WMAttachRightSide();

    WMAddEventHandler(
	IW_ALL_WINDOWS, EnterWindowMask, EventHandler, EnterWindowMask
    );
    WMAddEventHandler(
	IW_ALL_WINDOWS, LeaveWindowMask, EventHandler, LeaveWindowMask
    );
    WMAddEventHandler(
	IW_ALL_WINDOWS, PointerMotionMask, EventHandler, PointerMotionMask
    );
    WMEnableIWLEvent();
    WMDisplay();
    WMAppMainLoop();

    return 0;
}


static int quit(int* p_unused)
{
    long mask = EnterWindowMask | LeaveWindowMask | PointerMotionMask;

    WMCancelEventHandler(IW_ALL_WINDOWS, mask);
    WMCancelDisplayChange(IW_ALL_WINDOWS);
    exit(0);
    return 0;
}


static int EventHandler(int itype, XEvent* event)
{
    char sWid[4];
    float xyz[3];
    float dmin, dmax, dmean;
    int mode, iwid, iz, itime, mxyz[3];
    static int nxyz[3], iwave;
    static int use_overlay = 0;

    switch (itype) {
    case EnterWindowMask:
	iwid = IWXEvtToWid(event);
	(void) sprintf(sWid, "%d", iwid);
	IMAlPrt(0);
	IMOpen(1, sWid, "ro");
	IMRdHdr(1, nxyz, mxyz, &mode, &dmin, &dmax, &dmean);
	IWRtDisSec_ZWT(1, &iz, &iwave, &itime);
	if (IWRtNumOVColors() > 0) {
	    use_overlay = 1;
	    IWStartOverlay(1, IW_RED);
	} else {
	    use_overlay = 0;
	    IWBgnFastGr(1);
	}
	giGID = -1;
	giNPt = 0; 
	break;

    case LeaveWindowMask:
	IWGrRmProc(1);
	if (use_overlay) {
	    IWClearOverlay(1);
	} else {
	    IWEndFastGr(1);
	}
	IWDisplay(1);
	IMClose(1);
	giGID = -1;
	giNPt = 0;
	break;

    case PointerMotionMask:
	IWXEvtToData(1, iwave, event, xyz, &itime);
	if (xyz[0] >= 0 && xyz[0] < nxyz[0] &&
	    xyz[1] >= 0 && xyz[1] < nxyz[1]) {
	    if (giNPt < sizeof(gfPt) / sizeof(gfPt[0])) {
		if (giGID != -1) {
		    IWGrRmGrID(1, giGID);
		}
		gfPt[giNPt].x = xyz[0];
		gfPt[giNPt].y = xyz[1];
		++giNPt;
		/*
		 * Since the new graphic is the same as the old one with an
		 * addtional point, do not need to clear.  In the general
		 * case, however, the application would call IWClearOverlay
		 * (if using the overlay planes) or IWDisFastGr (if not) at
		 * this point.
		 */
		giGID = IWGrAddLns2D(1, gfPt, giNPt, 1, 0, IW_RED);
		IWGrDisGrID(1, giGID);
	    }
	}
	break;

    default:
	WMPostInfo("Unexpected event detected", 5);
	break;
    }

    return 0;
}

