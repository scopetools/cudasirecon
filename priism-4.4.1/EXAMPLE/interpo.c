/*
 * This application is essentially the same as Resample2D in Priism: it
 * interpolates the input data to a grid that may have been translated,
 * rotated, stretched, and skewed in two dimensions.
 *
 * The interaction with the IP framework is unusual because the output may
 * have a different x and y dimensions than the selected input region.  To
 * handle the different dimensions, a callback, ComputeSize, is registered
 * with IPSetComputeSize; that callback informs the framework what the
 * output size is.  Another callback, FixOutputHeader, is registered
 * with IPAddOutput.  This callback is called after the framework has opened
 * the output file and set its header; FixOutputHeader then adjusts the
 * the geometry description (pixel size, tilt angles, origin) in the header
 * to account for the resampling.  FixOutputHeader also has the side effect
 * of calculating the coordinate transformation matrix used to perform
 * the resampling; it would be natural to calculate that matrix in an INIT_PROC
 * routine but since some of the same calculations are necessary to update
 * the header, the matrix calculation is done in FixOutputHeader.
 *
 * Does not handle complex-valued data.
 */

#include "ip.h"
#include "WMInclude.h"
#include "IWInclude.h"
#include <stdio.h>      /* sprintf */
#include <stdlib.h>     /* free, malloc */
#include <string.h>     /* strcmp, strlen */


static void Verify(int unused);
static void ComputeSize(int out_stream, int unused, int sizes[5]);
static void FixOutputHeader(int unused);
static void NewArea(int unused);
static void NewInputFile(int unused);
static void ProcData(int nxyz[3], float dmmm[3], int unused);
static void PromptForParameters(int unused);
static void ReadCmdLine(int unused);
static void DisplayUsage(int unused);
static void WriteCmdLine(int unused);
static int NewMag(int* p_unused);
static void SetMenu(int unused);


static float gfRotX = 0.0f;
static float gfTranslation[2] = { 0.0f, 0.0f };
static float gfAngle = 90.0f;
static float gfMag[2] = { 1.0f, 1.0f };
static float gfRotCen[2] = { 0.0f, 0.0f };
static int giNewDim[2] = { 0, 0 }; 
static int giKeepCellDim = 0;
static float matrix[4];
static int giAutoSize = 1;


int
IPAppSpecifics(void)
{
    IPAddInput(1, " In:", True, False, NewInputFile, 0);
    IPAddOutput(2, "Out:", True, False, FixOutputHeader, 0);
    /*
     * Since complex values are not handled, do not allow the user to choose
     * complex output.
     */
    IPAllowMode(2, IW_COMPLEX, 0);
    IPAllowMode(2, IW_COMPLEX_SHORT, 0);
    /*
     * The processed output can have a different size than the input so notify
     * the IP framework.
     */
    IPSetComputeSize(2, ComputeSize, 0);
    IPSetCustRoutine(NEW_REGION, NewArea, 0);
    IPSetCustRoutine(CHECK_PARAMS, Verify, 0);
    IPSetCustRoutine(PROC_FUNC, ProcData, 0);
    IPSetMenuTitle("2D Interpolation");
    IPSetMenuLoc(450, 400);
    IPSetCustRoutine(CUS_MENU, SetMenu, 0);
    IPSetCustRoutine(CMD_LINE_FUNC, ReadCmdLine, 0);
    IPSetCustRoutine(USAGE_FUNC, DisplayUsage, 0);
    IPSetCustRoutine(PROMPT_FUNC, PromptForParameters, 0);
    IPSetCustRoutine(WRITE_CMD, WriteCmdLine, 0);

    return 0;
}


static
void
Verify(int unused)
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
    if (giNewDim[0] < 1 || giNewDim[1] < 1) {
	IPDisplayMessage("Output dimensions must be positive", IP_LOG_ERROR);
	IPDesist();
    }
}


static
void
ComputeSize(int out_stream, int unused, int sizes[5])
{
    sizes[0] = giNewDim[0];
    sizes[1] = giNewDim[1];
}


static
void
FixOutputHeader(int unused)
{
    int nxyz[3], mxyz[3], nxyzst[3];
    int i, j;
    float cell[6];
    float deltai[3], deltao[3];
    float new_orig[3];
    float tilts[3], oldrot[3][3], approt[3][3], newrot[3][3];

    IMRtSiz(1, nxyz, mxyz, nxyzst);
    IMRtCel(1, cell);
    deltai[0] = cell[0] / mxyz[0];
    deltai[1] = cell[1] / mxyz[1];
    deltai[2] = cell[2] / mxyz[2];
    deltao[0] = cell[0] / mxyz[0] / gfMag[0]; 
    deltao[1] = cell[1] / mxyz[1] / gfMag[1]; 
    deltao[2] = deltai[2];
    if (! giKeepCellDim) {
	int bounds[4][3], iwave[IW_MAX_WAVE];

	IPGetRegionInfo(bounds, iwave);
	cell[0] = deltao[0];
	cell[1] = deltao[1];
	cell[2] = deltao[2] * bounds[2][2];
	if (mxyz[0] != 1 || mxyz[1] != 1 || mxyz[2] != 1) {
	    cell[0] *= giNewDim[0];
	    cell[1] *= giNewDim[1];
	    cell[2] *= (bounds[2][1] - bounds[2][0]) /
		((bounds[2][2] != 0) ? bounds[2][2] : 1) + 1;
	}
    }
    cell[5] = gfAngle; 

    if (! IPIsAppending(2)) {
	IMRtSiz(2, nxyz, mxyz, nxyzst);
	/*
	 * Because the data has been resampled, rotated, and or skewed relative
	 * to the original, it is not clear how to interpret the starting
	 * indices.  Simply setting them to zero.
	 */
	nxyzst[0] = 0; nxyzst[1] = 0; nxyzst[2] = 0;
	IMAlSiz(2, nxyz, nxyzst);
	IMAlCel(2, cell);

	IMRtTlt(1, tilts);
	IWMatrixFromAngles(tilts, oldrot);
	tilts[0] = 0.0f;
	tilts[1] = 0.0f;
	tilts[2] = -gfRotX;
	IWMatrixFromAngles(tilts, approt);
	for (j = 0; j < 3; ++j) {
	    for (i = 0; i < 3; ++i) {
		newrot[j][i] = approt[0][i] * oldrot[j][0] +
		    approt[1][i] * oldrot[j][1] + approt[2][i] * oldrot[j][2];
	    }
	}
	IWAnglesFromMatrix(tilts, newrot);
	IMAlTlt(2, tilts);

	if (gfAngle == 90.0f) {
	    int nxyzt[4][3], iwave[IW_MAX_WAVE];
	    float orig[3];

	    IPGetRegionInfo(nxyzt, iwave);
	    IMRtOrig(1, &orig[0], &orig[1], &orig[2]);
	    for (j = 0; j < 2; ++j) {
		new_orig[j] = -(gfTranslation[j] + 0.5 * giNewDim[j]) *
		    deltao[j];
		for (i = 0; i < 2; ++i) {
		    new_orig[j] += approt[j][i] *
			((gfRotCen[i] + nxyzt[i][0]) * deltai[i] + orig[i]);
		}
		new_orig[j] += approt[j][2] *
		    ((0.5 * ((nxyzt[2][1] - nxyzt[2][0]) / nxyzt[2][2] + 1) +
		      nxyzt[2][0]) * deltai[2] + orig[2]);
	    }
	    new_orig[2] = orig[2] + nxyzt[2][0] * deltai[2];
	} else {
	    /*
	     * The library coordinate transformations (IWDataToDVPos
	     * for instance) don't take into account the angles between
	     * the axes, so it doesn't particularly help to specify how
	     * to map the original coordinates into the new ones.
	     */
	    new_orig[0] = 0.0f;
	    new_orig[1] = 0.0f;
	    new_orig[2] = 0.0f;
	}
	IMAlOrig(2, new_orig[0], new_orig[1], new_orig[2]);
    }

    /*
     * Assumes that the Fortran compiler mangles names by appending an
     * underscore.
     */
    gettransmatrix_(
	&deltai[0],
	&deltai[1],
	&gfRotX,
	&gfAngle,
	&deltao[0],
	&deltao[1],
	matrix
    );
}


static
void
NewArea(int unused)
{
    int iwave[IW_MAX_WAVE];
    int nxyzt[4][3], nx, ny;

    IPGetRegionInfo(nxyzt, iwave);
    nx = nxyzt[0][1] - nxyzt[0][0] + 1;
    ny = nxyzt[1][1] - nxyzt[1][0] + 1;
    giNewDim[0] = nx * gfMag[0];
    giNewDim[1] = ny * gfMag[1];
    gfRotCen[0] = nxyzt[0][0] + nx * 0.5f;
    gfRotCen[1] = nxyzt[1][0] + ny * 0.5f;
    if (IPHasGUI()) {
	WMUpdateGroup(4);
    }
}


static
void
NewInputFile(int unused)
{
    float cell[6];

    IMRtCel(1, cell);
    gfAngle = cell[5];
    NewArea(0);
}


static
void
ProcData(int nxyz[3], float dmmm[3], int unused)
{
    int nxyzt[4][3], nx, ny;
    int iwave[IW_MAX_WAVE];
    float* buf = IPGetDataPtr(1);
    float* out = IPGetDataPtr(2);
    float cenx, ceny;

    IPGetRegionInfo(nxyzt, iwave);
    nx = nxyzt[0][1] - nxyzt[0][0] + 1;
    ny = nxyzt[1][1] - nxyzt[1][0] + 1;
    /*
     * IMInterp uses a different convention for the center so need
     * to compensate.
     */
    cenx = gfRotCen[0] - nxyzt[0][0] - 0.5f;
    ceny = gfRotCen[1] - nxyzt[1][0] - 0.5f;
    IMInterp(
	buf,
	out,
	nx,
	ny,
	giNewDim[0],
	giNewDim[1],
	matrix,
	cenx,
	ceny,
	gfTranslation[0] - 0.5f,
	gfTranslation[1] - 0.5f,
	1.0f
    );
    IMCalcDen(
	out,
	giNewDim[0],
	giNewDim[1],
	1,
	giNewDim[0],
	1,
	giNewDim[1],
	&dmmm[0],
	&dmmm[1],
	&dmmm[2]
    );
    dmmm[2] *= giNewDim[0] * giNewDim[1];
}


static
void
PromptForParameters(int unused)
{
    char yesno[8];

    (void) IPPromptReal("Z rotation", 1, 1, 1, &gfRotX);
    (void) IPPromptReal("X/Y shift", 2, 2, 2, gfTranslation);
    (void) IPPromptReal("Angle between axes", 1, 1, 1, &gfAngle);
    (void) IPPromptReal("X/Y magnification", 2, 2, 2, gfMag);
    NewArea(0);
    (void) IPPromptReal("Center", 2, 2, 2, gfRotCen);
    (void) IPPromptInt("Output X/Y size", 2, 2, 2, giNewDim);
    (void) IPPromptString(
	"Keep cell dimensions (y/n; default n)= ", yesno, sizeof(yesno)
    );
    if (yesno[0] == 'y' || yesno[0] == 'Y') {
	giKeepCellDim = 1;
    } else {
	giKeepCellDim = 0;
    }
}


static
void
ReadCmdLine(int unused)
{
    int i = 1;
    int size_set = 0;
    int center_set = 0;
    char* arg;

    giKeepCellDim = 0;

    while ((arg = IPGetArg(i)) != 0) {
	char leader[] = "Unrecognized option";
	char* msg;
	int ncnvt;

	++i;
	ncnvt = IPParseRealArg(arg, "rotation", 1, &gfRotX);
	if (ncnvt != 0) {
	    if (ncnvt != 1) {
		IPDisplayMessage(
		    "-rotation requires one real-valued argument", IP_LOG_ERROR
		);
		IPDesist();
	    }
	    continue;
	}
	ncnvt = IPParseRealArg(arg, "shift", 2, gfTranslation);
	if (ncnvt != 0) {
	    if (ncnvt != 2) {
		IPDisplayMessage(
		    "-shift requires two real-valued arguments", IP_LOG_ERROR
		);
		IPDesist();
	    }
	    continue;
	}
	ncnvt = IPParseRealArg(arg, "axes_angle", 1, &gfAngle);
	if (ncnvt != 0) {
	    if (ncnvt != 1) {
		IPDisplayMessage(
		    "-axes_angle requires one real-valued argument",
		    IP_LOG_ERROR
		);
		IPDesist();
	    }
	    continue;
	}
	ncnvt = IPParseRealArg(arg, "mag", 2, gfMag);
	if (ncnvt != 0) {
	    if (ncnvt != 2) {
		IPDisplayMessage(
		    "-mag requires two real-valued arguments", IP_LOG_ERROR
		);
		IPDesist();
	    }
	    continue;
	}
	ncnvt = IPParseRealArg(arg, "center", 2, gfRotCen);
	if (ncnvt != 0) {
	    if (ncnvt == 2) {
		center_set = 1;
	    } else {
		IPDisplayMessage(
		    "-center requires two real-valued arguments", IP_LOG_ERROR
		);
	    }
	    continue;
	}
	ncnvt = IPParseIntArg(arg, "size", 2, giNewDim);
	if (ncnvt != 0) {
	    if (ncnvt == 2) {
		size_set = 1;
	    } else {
		IPDisplayMessage(
		    "-size requires two integer arguments", IP_LOG_ERROR
		);
		IPDesist();
	    }
	    continue;
	}
	if (strcmp(arg, "-same_cell") == 0) {
	    giKeepCellDim = 1;
	    continue;
	}
	msg = (char*) malloc(sizeof(leader) + 1 + strlen(arg));
	if (msg != 0) {
	    (void) sprintf(msg, "%s %s", leader, arg);
	    IPDisplayMessage(msg, IP_LOG_ERROR);
	    free(msg);
	} else {
	    IPDisplayMessage(leader, IP_LOG_ERROR);
	}
	IPDesist();
    }

    if (! size_set || ! center_set) {
	int nxyzt[4][3], iwave[IW_MAX_WAVE];
	int nx, ny;

	IPGetRegionInfo(nxyzt, iwave);
	nx = nxyzt[0][1] - nxyzt[0][0] + 1;
	ny = nxyzt[1][1] - nxyzt[1][0] + 1;
	if (! size_set) {
	    giNewDim[0] = nx * gfMag[0];
	    giNewDim[1] = ny * gfMag[1];
	}
	if (! center_set) {
	    gfRotCen[0] = nxyzt[0][0] + nx * 0.5f;
	    gfRotCen[1] = nxyzt[1][0] + ny * 0.5f;
	}
    }
}


static
void
DisplayUsage(int unused)
{
    IPDisplayMessage("Application-specific arguments:", IP_LOG_ERROR);
    IPDisplayMessage(
	"     [-rotation=angle_deg] [-shift=x_shift:y_shift] \\", IP_LOG_ERROR
    );
    IPDisplayMessage(
	"     [-axes_angle=angle_deg] [-mag=x_mag:y_mag] \\", IP_LOG_ERROR
    );
    IPDisplayMessage(
	"     [-center=x_center:ycenter] [-size=nx_out:ny_out] [-same_cell]",
	IP_LOG_ERROR
    );
}


static
void
WriteCmdLine(int unused)
{
    char arg[60];

    (void) sprintf(arg, "-rotation=%.3f", gfRotX);
    IPAppendCommand(arg);
    (void) sprintf(
	arg, "-shift=%.1f:%.1f", gfTranslation[0], gfTranslation[1]
    );
    IPAppendCommand(arg);
    (void) sprintf(arg, "-axes_angle=%.3f", gfAngle);
    IPAppendCommand(arg);
    (void) sprintf(arg, "-mag=%.3f:%.3f", gfMag[0], gfMag[1]);
    IPAppendCommand(arg);
    (void) sprintf(arg, "-center=%.1f:%.1f", gfRotCen[0], gfRotCen[1]);
    IPAppendCommand(arg);
    (void) sprintf(arg, "-size=%d:%d", giNewDim[0], giNewDim[1]);
    IPAppendCommand(arg);
    if (giKeepCellDim) {
	IPAppendCommand("-same_cell");
    }
}


static
int
NewMag(int* p_unused)
{
    if (giAutoSize) {
	NewArea(0);
    }
    return 0;
}


static
void
SetMenu(int unused)
{
    static float maxrotx = 180.0f;
    static float minrotx = -180.0f;
    static float fzero = 0.0f;

    WMAddInfoButton("         Z rotation angle:", "interpo xrot");
    WMAddFloatField(
	&gfRotX, 1, 8, 16, &minrotx, &maxrotx, &three, NULL, NULL, 0, 0
    );
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("          X,Y translation:", "interpo translation");
    WMAddFloatField(
	gfTranslation, 2, 9, 32, &fzero, &fzero, &one, NULL, NULL, 0, 0
    );
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("       Angle between axes:", "interpo Angle");
    WMAddFloatField(
	&gfAngle, 1, 8, 16, &fzero, &fzero, &three, NULL, NULL, 0, 4
    );
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("X, Y magnification factor:", "interpo mag");
    WMAddFloatField(
	gfMag, 2, 11, 32, &fzero, &fzero, &three, NewMag, NULL, 0, 0
    );
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("    Transformation center:", "interpo center");
    WMAddFloatField(
	gfRotCen, 2, 10, 32, &fzero, &fzero, &one, NULL, NULL, 0, 4
    );
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("         Output X, Y size:", "interpo size");
    WMAddIntField(giNewDim, 2, 10, 15, &one, &one, &one, NULL, NULL, 0, 4);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("Same Cell dim:", "interpo cell_dim");
    WMAddToggleButton("", &giKeepCellDim, NULL, NULL, 0, 4);
    WMAddInfoButton("Autosize:", "interpo autosize");
    WMAddToggleButton("", &giAutoSize, NULL, NULL, 0, 4);
    WMNewRow();
}
