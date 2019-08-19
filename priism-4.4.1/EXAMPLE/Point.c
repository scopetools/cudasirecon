/*
 * This application is essentially Priism's Point Values, but the extra
 * complexities of handling the cross-section profiles have been removed.
 * It monitors all image windows for events and the  IWLockSec_ZWT and
 * IWUnlockSec_ZWT to gain zero-copy access to the image data.  Since it
 * is only accessing a single pixel value at a time, it is arguably better to
 * use IMRdPas to retrieve the data value since that avoids the complexities
 * of dealing with all the possible pixel formats (see FindIntensities) and if
 * the data had been swapped out, it avoids swapping in the entire section's
 * worth of image data.
 */

#include "IWInclude.h"
#include "WMInclude.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#define POS_INTENS_GROUP (1 << 8)
#define FIELD_WIDTH 21


static int EnterWindow(int* w, XEvent* event);
static void FindIntensities(
    const void* pData,
    int mode,
    int complexDisplayMode,
    const int ixyz[3],
    int nx,
    float* intensity
);
static int HandleDisplayChange(int unused, XEvent* event);
static int LeaveWindow(int* w, XEvent* event);
static int Motion(int* w, XEvent* event);
static int quit(int* pUnused);
static int quitWrapper(void);
static void ResetSizes(void);
static void UpdateDisplay(void);


static const int gInputStream = 1;
static int gWindowOpen = 0;

static char gWindowStr[8] = "     ";
static Widget gWindowField;

static char* gWaveLabels[IW_MAX_WAVE] = {" 0 ", " 1 ", " 2 ", " 3 ", " 4 "};
static Widget gWaveMenu;
static int gWaveShown = 0;

static char gCoordStr[4][32] = {
  "                              ",
  "                              ",
  "                              ",
  "                              "
};
static char gIntenStr[32] = {
  "                              "
};

static int gNXYZTW[5] = { 0, 0, 0, 0, 0 };
static int gMode = 0;

static float gCursorXYZ[3];
static int gCursorTime;


int
main(int argc, char* argv[])
{
    IMAlPrt(0);

    WMInit("Point Info");
    WMSetLoc(512, 512);
    WMAddInfoButton("Window", "Point WINDOW");
    gWindowField = WMAddText(gWindowStr, sizeof(gWindowStr) - 1, 1);
    WMAddInfoButton("Wave", "Point WAVE_NUM");
    gWaveMenu = WMAddOptionMenu(
	gWaveLabels, IW_MAX_WAVE, &gWaveShown, NULL, NULL, 0, 0
    );
    WMNewRow();
    WMSetOffset(0, 0, 8, 0);
    WMAddInfoButton("   X  ", "Point X_Y");
    WMAddText(gCoordStr[0], FIELD_WIDTH, POS_INTENS_GROUP);
    WMAttachRightSide();
    WMSetOffset(0, 0, 0, 0);
    WMNewRow();
    WMAddInfoButton("   Y  ", "Point X_Y");
    WMAddText(gCoordStr[1], FIELD_WIDTH, POS_INTENS_GROUP);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("   Z  ", "Point X_Y");
    WMAddText(gCoordStr[2], FIELD_WIDTH, POS_INTENS_GROUP);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("   T  ", "Point X_Y");
    WMAddText(gCoordStr[3], FIELD_WIDTH, POS_INTENS_GROUP);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("Inten.", "Point Intensity");
    WMAddText(gIntenStr, FIELD_WIDTH, POS_INTENS_GROUP);
    WMAttachRightSide();

    WMNewRow();
    WMAddInfoButton(" ? ", "Point");
    WMAddFuncButton("Quit", quit, NULL, 0, 0); 
    WMSetExitFunction(quitWrapper);
    WMAttachRightSide();

    WMDisplay();

    if (WMAddEventHandler(IW_ALL_WINDOWS, EnterWindowMask, EnterWindow, 0) !=
	IW_ERROR &&
	WMAddEventHandler(IW_ALL_WINDOWS, LeaveWindowMask, LeaveWindow, 0) !=
	IW_ERROR &&
	WMProcDisplayChange(IW_ALL_WINDOWS, HandleDisplayChange, 0) !=
	IW_ERROR) {
	WMEnableIWLEvent();
	WMAppMainLoop();
    } else {
	WMCancelEventHandler(
	    IW_ALL_WINDOWS, EnterWindowMask | LeaveWindowMask
	);
	WMConfirmError("Point Values is unable to attach to IVE session.");
    }

    return 1;
}


static
int
EnterWindow(int* w, XEvent* event)
{
    int window = IWIsIWLEvt(event);

    (void) LeaveWindow(NULL, NULL);

    if (window > 0 && IWAttachWin(gInputStream, window, "ro") != IW_ERROR) {
	int iw;

	gWindowOpen = window;

	(void) sprintf(gWindowStr, "%5d", window);
	WMUpdateField(gWindowField);

	WMAddEventHandler(gInputStream, PointerMotionMask, Motion, 0);

	ResetSizes();
	IWXEvtToDataExt(
	    gInputStream, gWaveShown, event, gCursorXYZ, &iw, &gCursorTime
	);
	if (iw != gWaveShown && iw >= 0 && iw < IW_MAX_WAVE && 
	    iw < gNXYZTW[4]) {
	    gWaveShown = iw;
	    WMUpdateField(gWaveMenu);
	}
	UpdateDisplay();
    }

    return 1;
}


static
void
FindIntensities(
    const void* pRaw,
    int mode,
    int complexDisplayMode,
    const int ixyz[3],
    int nx,
    float* intensity
)
{
    switch (mode) {
    case IW_BYTE:
	{
	    const unsigned char* pData = (const unsigned char*) pRaw;

	    *intensity = pData[ixyz[0] + nx * ixyz[1]];
	}
	break;

    case IW_SHORT:
	/* Fall through. */
    case IW_EMTOM:
	{
	    const short* pData = (const short*) pRaw;

	    *intensity = pData[ixyz[0] + nx * ixyz[1]];
	}
	break;

    case IW_FLOAT:
	{
	    const float* pData = (const float*) pRaw;

	    *intensity = pData[ixyz[0] + nx * ixyz[1]];
	}
	break;

    case IW_COMPLEX_SHORT:
	{
	    const short* pData = (const short*) pRaw;
	    int iRealPart;

	    iRealPart = 2 * (ixyz[0] + ixyz[1] * nx);
	    switch (complexDisplayMode) {
	    case IW_COMPLEX_AMPLITUDE:
		*intensity = sqrt(pData[iRealPart] * pData[iRealPart] +
				  pData[iRealPart + 1] * pData[iRealPart + 1]);
		break;

	    case IW_COMPLEX_REAL:
		*intensity = pData[iRealPart];
		break;

	    case IW_COMPLEX_IMAGINARY:
		*intensity = pData[iRealPart + 1];
		break;

	    case IW_COMPLEX_PHASE:
		if (pData[iRealPart] == 0.0f && pData[iRealPart + 1] == 0.0f) {
		    *intensity = 0.0f;
		} else {
		    *intensity = atan2(pData[iRealPart + 1], pData[iRealPart]);
		}
		break;
	    }
	}
	break;

    case IW_COMPLEX:
	{
	    const float* pData = (const float*) pRaw;
	    int iRealPart;

	    iRealPart = 2 * (ixyz[0] + ixyz[1] * nx);
	    switch (complexDisplayMode) {
	    case IW_COMPLEX_AMPLITUDE:
		*intensity = sqrt(pData[iRealPart] * pData[iRealPart] +
				  pData[iRealPart + 1] * pData[iRealPart + 1]);
		break;

	    case IW_COMPLEX_REAL:
		*intensity = pData[iRealPart];
		break;

	    case IW_COMPLEX_IMAGINARY:
		*intensity = pData[iRealPart + 1];
		break;

	    case IW_COMPLEX_PHASE:
		if (pData[iRealPart] == 0.0f && pData[iRealPart + 1] == 0.0f) {
		    *intensity = 0.0f;
		} else {
		    *intensity = atan2(pData[iRealPart + 1], pData[iRealPart]);
		}
		break;
	    }
	}
	break;

    case IW_USHORT:
	{
	    const unsigned short* pData = (const unsigned short*) pRaw;

	    *intensity = pData[ixyz[0] + nx * ixyz[1]];
	}
	break;

    case IW_LONG:
	{
	    const int* pData = (const int*) pRaw;

	    *intensity = pData[ixyz[0] + nx * ixyz[1]];
	}
	break;

    case IW_U4BIT:
	{
	    const unsigned char* pData = (const unsigned char*) pRaw;
	    int nx_round = (nx + 1) / 2;

	    if (ixyz[0] % 2 == 0) {
		*intensity =
		    pData[ixyz[0] / 2 + (size_t) nx_round * ixyz[1]] & 0xf;
	    } else {
		*intensity =
		    (pData[ixyz[0] / 2 + (size_t) nx_round * ixyz[1]] &
		     0xf0) >> 4;
	    }
	}
	break;
    }
}


static
int
HandleDisplayChange(int unused, XEvent* event)
{
    int affectedWindow = IWIsIWLEvt(event);

    if (affectedWindow == gWindowOpen) {
	switch (event->xclient.data.s[IW_CM_CHANGE_OR_SYN]) {
	case IW_WIN_KILLED:
	    (void) LeaveWindow(NULL, NULL);
	    break;

	default:
	    ResetSizes();
	    UpdateDisplay();
	    break;
	}
    }

    return 1;
}


static
int
LeaveWindow(int *w, XEvent *event)
{
    if (gWindowOpen) {
	int iCoord;

	WMCancelEventHandler(gInputStream, PointerMotionMask);

	(void) sprintf(gWindowStr, "%5s", " ");
	for (iCoord = 0; iCoord < 4; ++iCoord) {
	    (void) sprintf(gCoordStr[iCoord], "%*s", FIELD_WIDTH, " ");
	}
	(void) sprintf(gIntenStr, "%*s", FIELD_WIDTH, " ");
	WMUpdateField(gWindowField);
	WMUpdateGroup(POS_INTENS_GROUP);

	IMClose(gInputStream);
	gWindowOpen = 0;
    }
    return 1;
}


static
int
Motion(int* w, XEvent* event)
{
    int iw;

    IWXEvtToDataExt(
	gInputStream, gWaveShown, event, gCursorXYZ, &iw, &gCursorTime
    );
    if (iw != gWaveShown && iw >= 0 && iw < IW_MAX_WAVE && iw < gNXYZTW[4]) {
	gWaveShown = iw;
	WMUpdateField(gWaveMenu);
    }
    UpdateDisplay();

    return 1;
}


static
int
quit(int* pUnused)
{
    WMCancelDisplayChange(IW_ALL_WINDOWS);
    WMCancelEventHandler(IW_ALL_WINDOWS, EnterWindowMask | LeaveWindowMask);
    (void) LeaveWindow(NULL, NULL);
    exit(0);

    return 1;
}


static
int
quitWrapper(void)
{
    return quit(NULL);
}


static
void
ResetSizes(void)
{
    int mxyz[3], interleaving;
    int izwt[3];
    float dmin, dmax, dmean;
    int iWave;

    IMRdHdr(gInputStream, gNXYZTW, mxyz, &gMode, &dmin, &dmax, &dmean);
    IMRtZWT(
	gInputStream, &gNXYZTW[2], &gNXYZTW[4], &gNXYZTW[3], &interleaving
    );

    IWRtDisSec_ZWT(gInputStream, &izwt[0], &izwt[1], &izwt[2]);
    gCursorXYZ[2] = izwt[0];
    gCursorTime = izwt[2];
    if (IWRtColorMode(gInputStream) == IW_PSEUDO) {
	gWaveShown = izwt[1];
    } else {
	if (gWaveShown >= gNXYZTW[4] ||
	    ! IWIsWaveMap(gInputStream, gWaveShown)) {
	    gWaveShown = 0;
	    for (iWave = 0;
		 iWave < gNXYZTW[4] && iWave < IW_MAX_WAVE;
		 ++iWave) {
		if (IWIsWaveMap(gInputStream, iWave)) {
		    gWaveShown = iWave;
		    break;
		}
	    }
	}
    }
    WMUpdateField(gWaveMenu);
}


static
void
UpdateDisplay(void)
{
    int ixyz[3];

    ixyz[0] = floor(gCursorXYZ[0]);
    ixyz[1] = floor(gCursorXYZ[1]);
    ixyz[2] = floor(gCursorXYZ[2]);
    if (ixyz[0] >= 0 && ixyz[0] < gNXYZTW[0] &&
	ixyz[1] >= 0 && ixyz[1] < gNXYZTW[1] &&
	ixyz[2] >= 0 && ixyz[2] < gNXYZTW[2] &&
	gCursorTime >= 0 && gCursorTime < gNXYZTW[3]) {
	int iCoord;
	int iSection;
	float realXYZ[3];
	float intensity;
	void* pData;

	IMRtSecNum(gInputStream, ixyz[2], gWaveShown, gCursorTime, &iSection);

	pData = IWLockSec_ZWT(gInputStream, ixyz[2], gWaveShown, gCursorTime);
	FindIntensities(
	    pData,
	    gMode,
	    IWRtComplexDis(gInputStream),
	    ixyz,
	    gNXYZTW[0],
	    &intensity
	);
	IWUnlockSec_ZWT(gInputStream, ixyz[2], gWaveShown, gCursorTime);

	IWDataToDVPos(
	    gInputStream, gCursorXYZ[0], gCursorXYZ[1], gCursorXYZ[2], realXYZ
	);
	for (iCoord = 0; iCoord < 3; ++iCoord) {
	    (void) sprintf(
		gCoordStr[iCoord],
		" %7.1f  %10.2f ",
		gCursorXYZ[iCoord],
		realXYZ[iCoord]);
	}
	(void) sprintf(
	    gCoordStr[3],
	    " %7.1f  %10.2f ",
	    (float) gCursorTime,
	    (float) gCursorTime
	);
	(void) sprintf(gIntenStr, "   %14.6g    ", intensity);
    } else {
	int iCoord;

	for (iCoord = 0; iCoord < 4; ++iCoord) {
	    (void) sprintf(gCoordStr[iCoord], "%*s", FIELD_WIDTH, " ");
	}
	(void) sprintf(gIntenStr, "%*s", FIELD_WIDTH, " ");
    }

    WMUpdateGroup(POS_INTENS_GROUP);
}
