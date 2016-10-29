/**************************** ProcWinEvent.c *********************************
   Description:
      A simple program to show how to implement data picking function in IVE
      by calling simple WML calls

   File Group: 
	ProcWinEvent.c 

   Author:   
        Hans Chen   1995 

   Rev 1.1:  reformat to IVE source file format      04/03/96
******************************************************************************/

#include "IWInclude.h"
#include "WMInclude.h"
#include <stdio.h>       /* sprintf */
#include <stdlib.h>      /* exit */
#include <string.h>      /* memset, strcpy */

      /* public function */

      /* local function */
static int quit(int* p_unused);
static int EventHandler(int itype, XEvent *event);

     /* public global variables */

     /* private global variable */
static char gcEventType[32];
static char giWid[4];
static int giStream = 1;
static int giWave = 0;
static IW_POINT gPtList[100];
static int giNumPt = 0;
static int giGID = -1, giCirGID = -1; 
static int giPolygon = 0;


int main(int argc, char* argv[])
{
    (void) memset(gcEventType, ' ', sizeof(gcEventType) - 1);
    gcEventType[sizeof(gcEventType) - 1] = '\0';
    (void) memset(giWid, ' ', sizeof(giWid) - 1);
    giWid[sizeof(giWid) - 1] = '\0';

    WMInit(NULL);
    WMSetLoc(500, 500);
    WMAddInfoButton("Window Event Demo", NULL);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("   Window ID", NULL);
    WMAddText(giWid, sizeof(giWid) - 1, 1);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("  Event info", NULL);
    WMAddText(gcEventType, sizeof(gcEventType) - 1, 1);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("Draw Polygon", NULL);
    WMAddToggleButton("", &giPolygon, NULL, NULL, 0, 0); 
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
    WMAddEventHandler(
	IW_ALL_WINDOWS, PointerMotionMask,EventHandler,PointerMotionMask
    );
    WMAddEventHandler(
	IW_ALL_WINDOWS, KeyPressMask, EventHandler, KeyPressMask
    );
    WMProcDisplayChange(IW_ALL_WINDOWS, EventHandler, -1);
    WMEnableIWLEvent();
    
    WMDisplay();
    WMAppMainLoop();

    return 0;
}


static int quit(int* p_unused)
{
    long mask = EnterWindowMask | LeaveWindowMask| ButtonPressMask|
	ButtonReleaseMask | PointerMotionMask;

    WMCancelEventHandler(IW_ALL_WINDOWS, mask);
    WMCancelDisplayChange(IW_ALL_WINDOWS);
    exit(0);
    return 0;
}

  
/*
 * itype : event type
 * event : pointer to X Event struct
 */
int EventHandler(int itype, XEvent *event)
{
    int iwid;
    int ixyz[3];
    float xyz[3];
    int itime;

    iwid = IWXEvtToWid(event);
    (void) sprintf(giWid, "%d", iwid);

    switch (itype) {
    case EnterWindowMask:
	IMAlPrt(0);
	IMOpen(giStream, giWid, "ro");
	giGID = -1;
	giCirGID = -1;
	(void) strcpy(gcEventType, "EnterWindowMask");
	break;

    case KeyPressMask:
	(void) sprintf(gcEventType, "keysym= %d", IWXEvtToKeySym(event));
	break;

    case LeaveWindowMask:
	giNumPt = 0;
	if (giGID != -1) {
	    IWGrRmGrID(giStream, giGID);
	    giGID = -1;
	}
	if (giCirGID != -1) {
	    IWGrRmGrID(giStream, giCirGID);
	    giCirGID = -1;
	}
	IWDisplay(giStream);
	IMClose(giStream);
	(void) strcpy(gcEventType, "LeaveWindowMask");
	break;

    case ButtonPressMask:
	{
	    int boxSize[] = {1, 1};
	    int icolor = IW_RED;
	    unsigned char pixel;
	    float val;
	    float scale[4];

	    IWXEvtToData(giStream, giWave, event, xyz, &itime);
	    ixyz[0] = xyz[0];
	    ixyz[1] = xyz[1];
	    ixyz[2] = xyz[2];
	    IMPosnZWT(giStream, ixyz[2], giWave, itime);
	    IMRdPas(giStream, &val, 1, 1, ixyz[0], ixyz[0], ixyz[1], ixyz[1]);
	    IWRtScl(giStream, giWave, scale);
	    IWScaleImage(&val, boxSize, IW_FLOAT, scale, &pixel);
	    if (giNumPt == sizeof(gPtList) / sizeof(gPtList[0])) {
		giNumPt = 0;
	    }
	    gPtList[giNumPt].x = ixyz[0];
	    gPtList[giNumPt].y = ixyz[1];
	    gPtList[giNumPt].z = ixyz[2];
	    giNumPt++;
	    if (giGID != -1) {
		IWGrRmGrID(giStream, giGID);
	    }
	    if (giCirGID != -1) {
		IWGrRmGrID(giStream, giCirGID);
	    }
	    giCirGID=IWGrAddCir2D(
		giStream, ixyz[0], ixyz[1], 2.0f, 1, icolor
	    );
	    if (giNumPt >= 2) {
		if (! giPolygon) {
		    giGID = IWGrAddLns3D(
			giStream, gPtList, giNumPt, 1, 1, icolor
		    );
		} else {
		    giGID = IWGrAddPoly3D(
			giStream, gPtList, giNumPt, 1, icolor
		    );
		}
	    } else {
		giGID = -1;
	    }
	    IWDisplay(giStream);
	    (void) sprintf(gcEventType, "Img:%g Pix:%d", val, (int)pixel);
	}
	break;

    case ButtonReleaseMask:
	(void) strcpy(gcEventType, "ButtonReleaseMask");
	break;

    case PointerMotionMask:
	IWXEvtToData(giStream, giWave, event, xyz, &itime);
	ixyz[0] = xyz[0];
	ixyz[1] = xyz[1];
	ixyz[2] = xyz[2];
	(void) sprintf(gcEventType, "X:%3d Y:%3d", ixyz[0], ixyz[1]);
	break;

    case -1:
	(void) sprintf(gcEventType, "Section #: %d\n", IWRtDisSec(giStream));
	break;

    default:
	(void) strcpy(gcEventType, "Unknown event");
    }

    WMUpdateGroup(1);

    return 0;
}


