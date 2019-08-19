/*
 * This is a simple program to show how to use WMAddScrolledText widget
 * to display more than one line of text  
 *
 * Hans Chen 
 */

#include "WMInclude.h"
#include <stdlib.h>  /* exit */
#include <string.h>  /* strcpy, memset */


char string[20];


int
quit()
{
   exit(0);
   return 0;
}


int
enter_text(char* str)
{
    (void) strcpy(string, str);
    WMUpdateGroup(1);
    (void) strcpy(str, "");
    WMUpdateGroup(2);
    return 0;
}


int
main(int argc, char* argv[])
{
    int nvisible_line = 3;
    int igroup = 1;
    int nsave = 5;
    char temp[sizeof(string)];

    (void) memset(string, ' ', sizeof(string) - 1);
    string[sizeof(string) - 1] = '\0';
    (void) memset(temp, ' ', sizeof(temp) - 1);
    temp[sizeof(temp) - 1] = '\0';

    WMInit(NULL);
    WMAddInfoButton("  Scroll / Text Widget  ", NULL);
    WMNewRow();
    WMAddScrolledText(string, &nsave, nvisible_line, igroup);
    WMNewRow();
    WMAddInfoButton("Text:", NULL);
    WMAddText(string, 15, igroup);
    WMAttachRightSide();
    WMNewRow();
    WMAddInfoButton("MSG:",NULL);
    WMAddCharField(temp, 20, 15, enter_text, temp, 0, 2);
    WMAttachRightSide();
    WMNewRow();
    WMAddFuncButton(" Exit ", quit, NULL, 0, 0);
    WMAttachRightSide();
    WMDisplay();
    WMAppMainLoop();

    return 0;
}
