/*
 * Simple program to show WMAddArrowButton function in WM.a
 * Hans 06/09/96
 */
#include "WMInclude.h"
#include <unistd.h>  /* sleep */
#include <stdio.h>   /* sprintf */
#include <stdlib.h>  /* exit */
#include <string.h>  /* memset */


char status[20];


int
callback(int *idir)
{
    (void) sprintf(status,"dir = %d but= %d", *idir, WMGetLastButNum());
    WMUpdateGroup(1);
    if (*idir == 1) {
	sleep(2);
	WMUndisplaySubMenu(0);
	sleep(5);
	WMDisplaySubMenu(0);
    }
    return 0;
}


int
quit()
{
    exit(0);
    return 0;
}
  

int
main(int argc, char* argv[])
{
    int one=1, two=2, three=3, four=4; 

    (void) memset(status, ' ', sizeof(status) - 1);
    status[sizeof(status) - 1] = '\0';

    WMInit(NULL);
    WMSetLoc(500, 500);
    WMAddInfoButton("Arrow Widget Demo", NULL);
    WMAttachRightSide();
    WMNewRow();
    WMAddArrowButton(WM_UP_ARROW, callback, &one);
    WMAddArrowButton(WM_DOWN_ARROW, callback, &two);
    WMAddArrowButton(WM_LEFT_ARROW, callback, &three);
    WMAddArrowButton(WM_RIGHT_ARROW, callback, &four);
    WMAddInfoButton("Status:", NULL);
    WMAddText(status, sizeof(status), 1);
    WMNewRow();
    WMAddFuncButton("Exit", quit, NULL, 0, 0);
    WMAttachRightSide();
    WMDisplay();
    WMAppMainLoop();

    return 0;
}

