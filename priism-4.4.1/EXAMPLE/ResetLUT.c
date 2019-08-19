/* A simple program to reset IWL pseudo color lut table  */


#include "IWInclude.h"
#include <stdio.h>    /* fprintf, printf, stderr */
#include <stdlib.h>   /* free, malloc */


int 
main(int argc, char* argv[])
{
    int imin, imax, num_values, csiz;
    int* ctab;
    int* ired;
    int* igreen;
    int* iblue;
    float scale;
    int i;

    if (IWAttach() == IW_ERROR) {
	(void) fprintf(
	    stderr,
	    "IWL not running or unable to attach to shared memory; quitting\n"
	);
	return 1;
    }
    
    num_values = IWRtLutSize();
    IWRtSclMinMax(&imin, &imax);
    (void) printf(" imin max = %d %d %d\n", imin, imax, num_values);
    csiz = imin + num_values + IW_NUM_COLORS;
    ctab = (int*) malloc(3 * sizeof(int) * csiz);
    if (ctab == 0) {
	(void) fprintf(stderr, "Unable to allocate space for color table\n");
	return 1;
    }
    ired = ctab;
    igreen = ctab + csiz;
    iblue = ctab + csiz * 2;

    IWRtLut(ired, igreen, iblue);
    scale = 255.0f / (imax - imin);
    for (i = imin; i <= imax; ++i) {
	int cindex = (i - imin) * scale;

	ired[i] = cindex;
	igreen[i] = cindex;
	iblue[i] = cindex;
    }
    IWAlLut(ired, igreen, iblue);

    free(ctab);

    return 0;
}
