/* a simple program allows user to pick a new image display center */

#include "IWInclude.h"
#include "WMInclude.h"
#include <stdio.h>      /* fprintf, sscanf */
#include <stdlib.h>     /* exit */
#include <string.h>     /* strlen */


static const int winStream = 1;
int giWave[IW_MAX_WAVE];


int NewCenter(int argu, XEvent* event);


int
main(int argc, char* argv[])
{
    static char msgPart1[] = "Use left mouse button to";
    static char msgPart2[] = "choose new display center";

    int iwid;
    int nxyz[3], mxyz[3], mode;
    float dmin, dmax, dmean;
    int nz, nw, nt, iflag;
    int window_x, window_y, window_width, window_height;
    int screen_width, screen_height;
    int iw;

    if (argc < 2 ||
	sscanf(argv[1], "%d", &iwid) != 1) {
	(void) fprintf(stderr, "usage: CenterImg windowNumber\n"); 
	exit(1);
    }
    if (! IWAttachWin(winStream, iwid, "ro")) {
	(void) fprintf(stderr, "CenterImg: failed to open window %d\n", iwid); 
	exit(1);
    }

    /*
     * Suppresses the printing of the image header information.
     */
    IMAlPrt(0);

    IMRdHdr(winStream, nxyz, mxyz, &mode, &dmin, &dmax, &dmean);
    IMRtZWT(winStream, &nz, &nw, &nt, &iflag);
    if (nw == 0) {
	nw = 1;
    }

    for (iw = 0; iw < nw; ++iw) {
	giWave[iw] = (IWIsWaveMap(winStream, iw)) ? 1 : 0;
    }
    for (; iw < IW_MAX_WAVE; ++iw) {
	giWave[iw] = 0;
    }

    WMInit(NULL);

    /*
     * Position the dialog next to the image window.  The coordinates
     * returned by IWRtDisArea are of the lower left corner and are relative
     * to the lower left corner of the screen.  Those need to be converted
     * into the coordinates expected by X (upper right corner relative to
     * the upper right of the screen). 
     *
     * WMGetScreenSize is a new call in 4.0.7.  In previous versions, most
     * applications would just hardwire a screen size.
     */
    IWRtDisArea(
	winStream, &window_x, &window_y, &window_width, &window_height
    );
    WMGetScreenSize(&screen_width, &screen_height);
    WMSetLoc(
	window_x + window_width, screen_height - window_y - window_height
    );

    WMAddInfoButton(" Image Display Center ", NULL );
    WMAttachRightSide();
    WMNewRow();
    if (nw > 1) {
	WMAddInfoButton("Wave:", NULL);
	for (iw = 0; iw < nw; ++iw) {
	    WMAddToggleButton("", &giWave[iw], NULL, NULL, 0, 0);
	}
	WMNewRow();
    }
    WMAddText(msgPart1, strlen(msgPart1), 0);
    WMAttachRightSide();
    WMNewRow();
    WMAddText(msgPart2, strlen(msgPart2), 0);
    WMAttachRightSide();

    WMAddEventHandler(winStream, ButtonPressMask, NewCenter, NULL);
    WMEnableIWLEvent();
    WMDisplay();
    WMAppMainLoop(); 

    /*
     * WMAppMainLoop never returns, but this suppresses a compiler
     * warning.
     */
    return 1;  
}


int
NewCenter(int argu, XEvent* event)
{
    float zoom, offx, offy;
    int ixst, iyst, iwidth, iheight;
    int ix, iy, ioff[2];
    int ndx, ndy;

    int iw;

    ix = event->xbutton.x;
    iy = event->xbutton.y;
    IWRtWinPos(winStream, &ixst, &iyst, &iwidth, &iheight);
    IWRtMulDis(winStream, &ndx, &ndy);
    zoom = IWRtZoom(winStream);

    iwidth = iwidth / ndx;
    iheight = iheight / ndy; 
    ix = ix - iwidth*(ix/iwidth); 
    iy = iy - iheight*(iy/iheight); 
    offx = ix - iwidth / 2.0 ; 
    offy = iheight/2  - iy ; 
    offx = offx / zoom; 
    offy = offy / zoom; 
    for (iw = 0 ; iw < IW_MAX_WAVE; ++iw) {
	if (giWave[iw]) {
	    IWRtDisOffset(winStream, iw, ioff);
	    ioff[0] = ioff[0] + offx;
	    ioff[1] = ioff[1] + offy;
	    IWAlDisOffset(winStream, iw, ioff);
	}
    }

    /*
     * Forces the image window to do a full clear and then redisplay.
     */
    IWAlClrBkg(winStream, 1);
    IWDisplay(winStream);

    WMCancelEventHandler(winStream, ButtonPressMask);
    WMDisableIWLEvent();
    exit(0);

    return 0;
}
