/* Demonstrates how to create an image file from scratch. */

#include "IWInclude.h"
#include <float.h>     /* FLT_MAX */

int main(int argc, char* argv[])
{
    int nxyz[3] = { 128, 128, 20 };
    int mxyz[3] = { 1, 1, 1 };
    char label[80] = "label";
    int count;
    short lilarray[128*128];
    int i, j, k;
    float dmin, dmax, dmean, secmean;

    /*
     * Open a new file using the name specified on the command line or
     * a default value if one was not specified.
     */
    IMOpen(1, (argc >= 2) ? argv[1] : "/var/tmp/junk.dat", "new");
    /* Create header for file. */
    IMCrHdr(1, nxyz, mxyz, IW_SHORT, label, 1);
    IMAlCon(1, 0);
    dmin = FLT_MAX;
    dmax = -FLT_MAX;
    dmean = 0.0f;
    IMWrHdr(1, label, 0, dmin, dmax, dmean);
    for (i = 0; i < nxyz[2]; ++i) {
	secmean = 0.0f;
	for (j = 0, count = 0; j < nxyz[1]; ++j) {
	    for (k = 0; k < nxyz[0]; ++k, ++count) {
		lilarray[count] = (i * 10) + k;
		if (dmin > lilarray[count]) {
		    dmin = lilarray[count];
		}
		if (dmax < lilarray[count]) {
		    dmax = lilarray[count];
		}
		secmean += lilarray[count]; 
	    }
	}
	secmean /= nxyz[0] * nxyz[1];
	dmean += secmean;
	IMWrSec(1, lilarray);     /* write the sections one by one */
    }
    dmean = dmean / (nxyz[2]);
    IMWrHdr(1, label, 1, dmin, dmax, dmean);
    IMClose(1);

    return 0;
}

