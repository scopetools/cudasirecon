/*
 * A simple program to show how to use lut related IWL routines to reverse
 * IWL's pseudo color LUT 
 */


#include "IWInclude.h"
#include <stdio.h>     /* fprintf, printf, stderr */
#include <stdlib.h>    /* free, malloc */

int 
main(int argc, char** argv)
{
    int imin, imax, csiz;
    int* ctab;
    int* ired;
    int* igreen;
    int* iblue;
    int i;

    if (IWAttach() == IW_ERROR) {
	(void) fprintf(
	    stderr,
	    "IWL not running or unable to attach to shared memory; quitting\n"
	);
	return 1;
    }

    IWRtSclMinMax(&imin, &imax);
    csiz = imin + IWRtLutSize() + IW_NUM_COLORS;
    ctab = (int*) malloc(3 * sizeof(int) * csiz);
    if (ctab == 0) {
	(void) fprintf(stderr, "Unable to allocate space for color table\n");
	return 1;
    }
    ired = ctab;
    igreen = ctab + csiz;
    iblue = ctab + csiz * 2;

    IWRtLut(ired, igreen, iblue);
    /* Print the components of the graphic colors. */
    for (i = imax + 1; i <= imax + IW_NUM_COLORS; ++i) {
	(void) printf(" rgb= %d %d %d %d\n", i, ired[i], igreen[i], iblue[i]);
    }
    /* Reverse the color table. */
    for (i = imin; i <= imax + IW_NUM_COLORS; ++i) {
	ired[i] = 255 - ired[i];
	igreen[i] = 255 - igreen[i];
	iblue[i] = 255 - iblue[i];
    }
    IWAlLut(ired, igreen, iblue);

    free(ctab);

    return 0;
}
