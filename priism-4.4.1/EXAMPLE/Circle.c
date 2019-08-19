/*
 * A simple ploting program to show mouse motion tracking function in IW
 * and WM libraries.
 */
#include "IWInclude.h"
#include "WMInclude.h"
#include <stdio.h>      /* sprintf */


int quit(int* p_unused);
int EventHandler(int itype, XEvent* event);


int
main(int argc, char* argv[])
{
    static char sMessage[30] = " Move mouse into a window ";
    static char sMessage1[30] = " and try mouse button click";

    WMInit("Draw Circle Demo");
    WMSetLoc(500, 500);
    WMAddText(sMessage, 29, 1); 
    WMAttachRightSide();
    WMNewRow();
    WMAddText(sMessage1, 29, 1); 
    WMAttachRightSide();
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
	IW_ALL_WINDOWS, ButtonPressMask, EventHandler, ButtonPressMask
    );
    WMAddEventHandler(
	IW_ALL_WINDOWS, ButtonReleaseMask, EventHandler, ButtonReleaseMask
    );
    WMEnableIWLEvent();
    WMDisplay();
    WMAppMainLoop();

    return 0;
}


int
quit(int* p_unused)
{
    long mask;

    mask = EnterWindowMask | LeaveWindowMask | ButtonPressMask |
	ButtonReleaseMask;
    WMCancelEventHandler(IW_ALL_WINDOWS, mask);
    WMCancelDisplayChange(IW_ALL_WINDOWS);
    exit(0);
    return 0;
}


int EventHandler(int itype, XEvent* event)
{
    int iwid;
    char sWid[3];
    float xyz[3];
    float dmin, dmax, dmean;
    int mode, id, itime;
    static int nxyz[3], mxyz[3], icircle = -1;
    static int rad = 1;

    if (itype == EnterWindowMask) {
	iwid = IWXEvtToWid(event);
	(void) sprintf(sWid, "%d", iwid);
	IMAlPrt(0);
	IMOpen(1, sWid, "ro");
	IMRdHdr(1, nxyz, mxyz, &mode, &dmin, &dmax, &dmean);
	rad = 1;
    } else if (itype == LeaveWindowMask && icircle != -1) {
	if (WMConfirm(" Remove circle ?")) {
	    IWGrRmProc(1);
	    IWDisplay(1);
	}
	IMClose(1);
    } else if (itype == ButtonPressMask) {
	if (icircle != -1) {
           IWGrRmGrID(1, icircle);
	}
	icircle = IWGrAddCir2D(
	    1,
	    nxyz[0] / 2.0f,
	    nxyz[1] / 2.0,
	    nxyz[0] / 2.0f,
	    rad,
	    rad % 10
	);
	IWGrDisGrID(1, icircle);
	IWXEvtToData(1, 0, event, xyz, &itime);
	id = IWGrAddCir2D(1, xyz[0], xyz[1], rad, 2, 2);
	IWGrDisGrID(1, id);
	++rad;
    }
    return 0;
}

