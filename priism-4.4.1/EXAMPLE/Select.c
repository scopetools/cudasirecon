/*
 * Copies a subset of the input dataset to the output dataset.  If the
 * input is complex but the output is not, the imaginary part is dropped.
 * If the input is real but the output is not, an imaginary part of 0 is
 * inserted.
 */

#include "ip.h"
#include "WMInclude.h"
#include "IWInclude.h"
#include <math.h>      /* sqrt */


static void InitProc(int unused);
static void ProcData(int nxyz[3], float dmmm[3], int unused);
static void DoNothing(int unused);


int g_in_complex = 0;
int g_out_complex = 0;


int IPAppSpecifics(void)
{
    /*
     * Simple copy, so there's no need to have separate arrays for the input
     * and output.
     */
    IPSetIncoreProc();
    IPAddInput(one, " In:", True, False, NULL, 0);
    IPAddOutput(two, "Out:", True, False, NULL, 0);
    IPSetCustRoutine(INIT_PROC, InitProc, 0);
    IPSetCustRoutine(PROC_FUNC, ProcData, 0);
    IPSetMenuTitle("Region Select");
    IPSetMenuLoc(450, 400);
    /*
     * This application is unusual because it does not require any additional
     * parameters beyond what the framework handles.  Callbacks which
     * do nothing are provided for CMD_LINE_FUNC, PROMPT_FUNC, and WRITE_CMD
     * to enable the command-line and prompt modes and to enable the graphical
     * user interface components for writing out a command file.
     */
    IPSetCustRoutine(CMD_LINE_FUNC, DoNothing, 0);
    IPSetCustRoutine(PROMPT_FUNC, DoNothing, 0);
    IPSetCustRoutine(WRITE_CMD, DoNothing, 0);
    return 0;
}


/*
 * Need to determine whether the input and output are complex.
 */
static void InitProc(int unused)
{
    int   nxyz[3], mxyz[3];
    int   pixel_type;
    float mmm[3];

    IMRdHdr(one, nxyz, mxyz, &pixel_type, &mmm[0], &mmm[1], &mmm[2]);
    g_in_complex = pixel_type == IW_COMPLEX || pixel_type == IW_COMPLEX_SHORT;
    IMRdHdr(two, nxyz, mxyz, &pixel_type, &mmm[0], &mmm[1], &mmm[2]);
    g_out_complex = pixel_type == IW_COMPLEX || pixel_type == IW_COMPLEX_SHORT;
}


void ProcData(int nxyz[3], float dmmm[3], int unused)
{
    int nxy = nxyz[0] * nxyz[1];
    float* des = (float*) IPGetDataPtr(one);
    int i;

    if (g_in_complex == g_out_complex) {
	/* No conversion necessary so just update the statistics. */
	if (g_out_complex) {
	    nxy *= 2;
	    for (i = 0; i < nxy; i += 2) {
		float amplitude =
		    sqrt(des[i] * des[i] + des[i + 1] * des[i + 1]);

		if (dmmm[0] > amplitude) {
		    dmmm[0] = amplitude;
		}
		if (dmmm[1] < amplitude) {
		    dmmm[1] = amplitude;
		}
		dmmm[2] += amplitude;
	    }
	} else {
	    for (i = 0; i < nxy; ++i) {
		if (dmmm[0] > des[i]) {
		    dmmm[0] = des[i];
		}
		if (dmmm[1] < des[i]) {
		    dmmm[1] = des[i];
		}
		dmmm[2] += des[i];
	    }
	}
    } else {
	int j;

	if (g_in_complex) {
	    /* Retains only the real part. */
	    for (i = 0, j = 0; i < nxy; ++i, j += 2) {
		des[i] = des[j];
		if (dmmm[0] > des[i]) {
		    dmmm[0] = des[i];
		}
		if (dmmm[1] < des[i]) {
		    dmmm[1] = des[i];
		}
		dmmm[2] += des[i];
	    }
	} else {
	    /* Set the imaginary part to zero. */
	    for (i = nxy - 1, j = 2 * nxy - 2; i >= 0; --i, j -= 2) {
		des[j] = des[i];
		des[j + 1] = 0.0f;
		if (dmmm[0] > des[j]) {
		    dmmm[0] = des[j];
		}
		if (dmmm[1] < des[j]) {
		    dmmm[1] = des[j];
		}
		dmmm[2] += des[j];		
	    }
	}
    }
}


static void DoNothing(int unused)
{
}
