/*
 * Given a dataset, returns a dataset where points which are above a threshold
 * but which have one 4-connected neighbor below the threshold are set to 1
 * and all other points are set to zero.
 *
 * Does not handle complex-valued data.
 */

#include "ip.h"
#include "WMInclude.h"
#include "IWInclude.h"
#include <stdio.h>     /* sprintf, sscanf */
#include <stdlib.h>    /* free, malloc */
#include <string.h>    /* strlen */


static void CheckParams(int unused);
static void ProcData(int nxyz[3], float dmmm[3], int unused);
static void PromptForParameters(int unused);
static void ReadCmdLine(int unused);
static void DisplayUsage(int unused);
static void WriteCmdLine(int unused);
static void SetMenu(int unused);
static int GetString(int* p_unused);


static float gfThreshold = 0.0f;
static char gsThreshold[24] = "0.0";


int IPAppSpecifics(void)
{
    IPAddInput(1, "Open", True, False, NULL, 0);
    IPAddOutput(2, "Save", True, False, NULL, 0);
    /*
     * Since complex values are not handled, do not allow the user to choose
     * complex output.
     */
    IPAllowMode(2, IW_COMPLEX, 0);
    IPAllowMode(2, IW_COMPLEX_SHORT, 0);
    IPSetCustRoutine(CHECK_PARAMS, CheckParams, 0);
    IPSetCustRoutine(PROC_FUNC, ProcData, 0);
    IPSetMenuTitle("Find Border");
    IPSetMenuLoc(450, 400);
    IPSetCustRoutine(CUS_MENU, SetMenu, 0);
    IPSetCustRoutine(CMD_LINE_FUNC, ReadCmdLine, 0);
    IPSetCustRoutine(USAGE_FUNC, DisplayUsage, 0);
    IPSetCustRoutine(PROMPT_FUNC, PromptForParameters, 0);
    IPSetCustRoutine(WRITE_CMD, WriteCmdLine, 0);
    IPEnableResolutionSelection();
    return 0;
}


static void CheckParams(int unused)
{
    if (IPIsOpen(1)) {
	int mode;

	IMRtMode(1, &mode);
	if (mode == IW_COMPLEX || mode == IW_COMPLEX_SHORT) {
	    IPDisplayMessage(
		"Can not handle complex-valued input", IP_LOG_ERROR
	    );
	    IPDesist();
	}
    }
}


static void ProcData(int nxyz[3], float dmmm[3], int unused)
{
    float* src = IPGetDataPtr(1);
    float* des = IPGetDataPtr(2);
    int i, ix, j, iy, imy, ipy, npix;

    imy = 0;
    iy = nxyz[0];
    ipy = nxyz[0] * 2;
    npix = 0;
    for (i = 0; i < nxyz[0]; ++i) {
	des[i] = 0.0f;
    }
    for (j = 1; j < nxyz[1] - 1; ++j) {       
       des[iy] = 0.0f;
       for (i = 1; i < nxyz[0] - 1; ++i) {
	   ix = i + iy;
	   des[ix] = 0.0f;
	   if (src[ix] >= gfThreshold) {
	       if (src[ix - 1] < gfThreshold || src[ix + 1] < gfThreshold ||
		   src[i + ipy] < gfThreshold || src[i + imy] < gfThreshold) {
		   des[ix] = 1.0f;
		   ++npix;
	       }
	   }
       }
       iy = iy + nxyz[0];
       des[iy - 1] = 0.0f;
       imy = imy + nxyz[0];
       ipy = ipy + nxyz[0];
    }
    for (i = iy; i < iy + nxyz[0]; ++i) {
	des[i] = 0.0f;
    }

    dmmm[0] = 0.0f;
    dmmm[1] = (npix == 0) ? 0.0f : 1.0f;
    dmmm[2] += npix;
}


static void PromptForParameters(int unused)
{
    (void) IPPromptReal("Threshold", 1, 1, 1, &gfThreshold);
}


static void ReadCmdLine(int unused)
{
    int i = 1;
    char* arg;

    while ((arg = IPGetArg(i)) != 0) {
	int status;

	++i;
	status = IPParseRealArg(arg, "threshold", 1, &gfThreshold);
	if (status == 1) {
	    (void) sprintf(gsThreshold, "%g", gfThreshold);
	} else if (status == 0) {
	    char leader[] = "Unrecognized option";
	    char* msg = (char*) malloc(sizeof(leader) + 1 + strlen(arg));

	    if (msg != 0) {
		(void) sprintf(msg, "%s %s", leader, arg);
		IPDisplayMessage(msg, IP_LOG_ERROR);
		free(msg);
	    } else {
		IPDisplayMessage(leader, IP_LOG_ERROR);
	    }
	    IPDesist();
	} else {
	    IPDisplayMessage(
		"-threshold requires one real-valued argument", IP_LOG_ERROR
	    );
	    IPDesist();
	}
    }
}


static void DisplayUsage(int unused)
{
    IPDisplayMessage(
	"Application-specific arguments:\n     [-threshold=value]",
	IP_LOG_ERROR
    );
}


static void WriteCmdLine(int unused)
{
    char arg[40];

    (void) sprintf(arg, "-threshold=%g", gfThreshold);
    IPAppendCommand(arg);
}


static void SetMenu(int unused)
{
    WMAddInfoButton("Threshold:", "FindBorder Threshold");
    WMAddCharField(gsThreshold, 5, sizeof(gsThreshold), GetString, NULL, 0, 4);
    WMAttachRightSide();
    WMNewRow();
}


static int GetString(int* p_unused)
{
    float dummy;
   
    if (sscanf(gsThreshold, "%g", &dummy) == 1) {
	gfThreshold = dummy;
    } else {
	(void) sprintf(gsThreshold, "%g", gfThreshold);
	WMUpdateGroup(4);
    }
    return 0;
}
