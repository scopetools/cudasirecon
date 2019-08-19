/*
 * WM widget library provides a set of simple function calls for application
 * programs to create a Motif style GUI(Graphical Users Interface) without
 * learning the Motif or X window programming. The three major steps in
 * creating a "form" style GUI are following:
 * 
 * 1) WMInit(title): char *title; string to be shown up in title bar Initialize
 * the WM library and create a Motif form widget as the base of all the
 * subsequently added widget.
 * 
 * 2) Adding widgets to the GUI and manipulation their row-column positions by
 * calling various WM library calls:
 * 
 * Widget WMAddText(string ,maxlen ,group) char *string;  the actual text will
 * be printed on the GUI Add a output only text widget to the GUI. It does
 * not expect any type of the input from users.
 * 
 * Widget WMAddFloatField(fptr ,num ,dislen ,maxlen ,min ,max ,decimal ,callback
 * ,para ,update ,group); float *fptr;    an array of float variables int
 * num;       number of elements in the float array int dislen;    number of
 * characters displayed on the screen int maxlen;    the maximum length of
 * the char string representing the entire array of float numbers float min
 * ,max; range of the accepted value; if max == min no ranged is imposed int
 * decimal; number of decimal points displayed for each number
 * int(*callback)();  an integer function to be called when any value in the
 * array is changed int *para;  a pointer to a argument array for the
 * callback function int update;  flag to indicate whether the content in
 * "grouped field" should be updated after callback is executed int group;
 * flag to indical this field should be part of the "grouped field"
 * 
 * Widget WMAddIntField(iptr ,num ,dislen ,maxlen ,min ,max ,factor , callback
 * ,para ,update ,group); int *iptr;    an array of integer variables int
 * num;       number of elements in the float array int dislen;    number of
 * characters displayed on the screen int maxlen;    the maximum length of
 * the char string representing the entire array of float numbers float min
 * ,max; range of the accepted value; if max == min no ranged is imposed int
 * factor; an integer value which is the common factor amoung all the values
 * specified; this is used as a criteria to reject unvalid values
 * int(*callback)();  an integer function to be called when any value in the
 * array is changed int *para;  a pointer to a argument array for the
 * callback function int update;  flag to indicate whether the content in
 * "grouped field" should be updated after callback is executed int group;
 * flag to indical this field should be part of the "grouped field"
 * 
 * Widget WMAddCharField(string ,maxlen ,dislen ,callback , para ,update
 * ,group); char *string;  pointer to a array of chars for receiving input
 * string.; it is calling routine's responsibility to alloc enough space for
 * the char array. The size should be equal to maxlen  int dislen;    number
 * of characters displayed on the screen int maxlen;    the maximum length of
 * the char string representing the entire array of float numbers
 * int(*callback)();  an integer function to be called when any value in the
 * array is changed int *para;  a pointer to a argument array for the
 * callback function int update;  flag to indicate whether the content in
 * "grouped field" should be updated after callback is executed int group;
 * flag to indical this field should be part of the "grouped field"
 * 
 * Widget WMAddToggleButton(label ,ivar ,callback ,para ,update ,group); char
 * *label;  toggle button labeling string int *ivar;   an input int pointer
 * to store toggle button status. on = > *ivar = 1; off = > *ivar = 0;
 * int(*callback)();  an integer function to be called when any value in the
 * array is changed int *para;  a pointer to a argument array for the
 * callback function int update;  flag to indicate whether the content in
 * "grouped field" should be updated after callback is executed int group;
 * flag to indical this field should be part of the "grouped field"
 * 
 * Widget WMAddPulldown(label ,nlabel ,ichoice ,callback ,para ,i update
 * ,group); char **label;  a list of labeling strings for pulldown choices
 * int nlabel;   number of choice in pulldown menu int *ichoice;  an input
 * int pointer to store the last choice which is numbered from top to button
 * in the pulldown list staring from 0  int(*callback)();  an integer
 * function to be called when any value in the array is changed int *para;  a
 * pointer to a argument array for the callback function int update;  flag to
 * indicate whether the content in "grouped field" should be updated after
 * callback is executed int group;  flag to indical this field should be part
 * of the "grouped field"
 * 
 * Widget WMAddFuncButton(label ,callback ,para ,update ,group); char *label;
 * button labeling string int(*callback)();  an integer function to be called
 * when button is pushed int *para;  a pointer to a argument array for the
 * callback function int update;  flag to indicate whether the content in
 * "grouped field" should be updated after callback is executed int group;
 * flag to indical this field should be part of the "grouped field"
 * 
 * Widget WMAddInfoButton(label ,keyword); char label;  button label string char
 * *keyword;  a char string used as an index for global on line help facility
 * 
 * Widget WMNewRow(); The widgets are organized into rows. The first widget
 * created will at the upper-left corner of the menu and the following
 * widgets will be append to its right until WMNewRow is called. Then a new
 * row starts. It assume widgets on the same row have the same height.
 * 
 * int WMSetLROffset(left_offset ,right_offset) The default gap between widget
 * on the same row is 0. This can be changed to left_offset+right_offset. The
 * values input will remain in effect until next WMSetLROffset call is made.
 * 
 * 
 *  WMDisplay() This is the routine for displaying GUI on the screen.
 * 
 */
#include "WMInclude.h"
#include <stdio.h>   /* printf, sprintf */
#include <stdlib.h>  /* exit */
#include <string.h>  /* strcmp, strcpy */


Widget          w1, w2;
int             num[2] = {2, 3};
float           rnum[2] = {2.23, 3.33};
char           *label[] = {"choice 1", "choice 2", "choice 3"};
char            string[10];
char            rtitle[40] = " This is a test menu";
int             toggle = 1;
char            dir[30] = "/usr/people/hans";
char            pattern[10] = "*.c";
char            filename[60] = "chr342.d3d";
char            blabel[10] = "button1";
char           *func1[] = {"func XXXXXXXXXX 1 ", "func 2", "func 3", "func 4"};
char           *func2[] = {"FUNC 1 ", "FUNC 2", "FUNC 3", "FUNC 4", "FUNC 5", "FUNC 6"};
char            batchjob[50] = "/usr/people/clyborne/IVE/WIDGET/batch";
char            message[40] = " message";
int             wsub;
int             ifunc = 1;
int             ichoice = 0;
int             islide = 30;


int 
pfilename()
{
    (void) printf(" filename = %s\n", filename);
    return 0;
}


int 
leave()
{
    exit(0);
    return 0;
}


int 
pback()
{
    (void) printf(" ichoice = %d\n", ichoice);
/*
    WMUnpasteWidget(w1);
*/
    strcpy(rtitle, "XXXXX XXXXX XXXXX XXXXX");
    return 0;

}


int 
dragfunc(int *val)
{
    (void) printf(" slide = %d  %d\n", *val, islide);
    *val = (*val) + 1;
    return 0;
}


int 
tcallback()
{
    num[0] = num[0] + 1;
    ichoice = (ichoice + 1) % 3;
    (void) strcpy(string, "BBBBB");
    (void) strcpy(rtitle, "AAAAA BBBB CCCC DDDD");
    return 0;
}


int
lcallback()
{
    (void) printf(" in lcallback %d\n", ifunc);
    WMPasteWidget(w1);
    if (ifunc % 2 == 0) {
	WMUpdateFuncList(0, func2, NULL, 6, 0);
	(void) strcpy(message, func2[0]);
    }
    else {
	WMUpdateFuncList(0, func1, NULL, 4, 0);
	(void) strcpy(message, func1[0]);
    }
    return 0;
}


int 
ccallback()
{
    num[1]++;
    return 0;
}


int
icallback()
{
    if (toggle) {
	toggle = 0;
    } else {
	toggle = 1;
    }
    strcpy(string, "AAAA");
    rnum[0] = rnum[0] * 1.5;
    return 0;
}


int 
rcallback()
{
    char string[20];

    (void) sprintf(string, " rnum = %.2f %.2f\n", rnum[0], rnum[1]);
    (void) printf(" post info\n");
    WMPostInfo(string, 5);
    return 0;
}


int 
bcallback()
{
    if (strcmp(blabel, "button1") == 0) {
	strcpy(blabel, "button2");
    } else {
	strcpy(blabel, "button1");
    }
    WMDisplaySubMenu(wsub);
    return 0;
}


int 
unmapsubmenu()
{
    WMUndisplaySubMenu(wsub);
    return 0;
}


int
main(int argc, char* argv[])
{
    int             imax = 10, imin = 0, ifactor = 1;
    float           fmax = 100.0, fmin = 0.0;
    int             decimal = 2;
    int             ilist = 10;
    static int      drag = 0;
    char            tlabel[10] = "test menu";

    WMInit("IVEMenu testing program");
    WMSetLoc(512, 512);
    w1 = WMAddGetFile(
	"Filename", dir, pattern, 60, 20, filename, pfilename, NULL, 1, 1
    );
    WMNewRow();
    WMAddInfoButton("Title", " Keyword string");
    strcpy(string, "XXXXX");
    WMAddCharField(string, 10, 5, &ccallback, NULL, 1, 1);
    WMAddPulldown(label, 3, &ichoice, pback, NULL, 1, 1);
    WMAddText(tlabel, 11, 4);
    WMNewRow();
    WMAddIntField(
	num, 2, 10, 10, &imin, &imax, &ifactor, icallback, NULL, 1, 1
    );
    WMAddToggleButton("Toggle", &toggle, tcallback, NULL, 1, 1);
    WMNewRow();
    WMAddFloatField(
	rnum, 2, 20, 20, &fmin, &fmax, &decimal, rcallback, NULL, 1, 1
    );
    WMNewRow();
    WMAddHorSlider(
	"Slider", 200, 0, 200, &islide, tcallback, NULL, dragfunc, &drag, 1, 1
    );
    WMNewRow();
    WMAddFuncList(func1, NULL, 4, 0, 100, &ifunc, lcallback, NULL, 1);
    WMNewRow();
    WMAddFuncButton("Exit", leave, NULL, 1, 1);
    WMAddFuncButton(blabel, &bcallback, NULL, 1, 1);
    WMAddFuncButton("DoIt", &bcallback, NULL, 1, 1);
    WMAddDoButton(NULL,NULL,"NULL", batchjob);
    WMAddSaveButton("menu.save");
    WMAddScrolledText(message, &ilist, 5, 1);
    wsub = WMInitSubMenu("1st submenu ");
    WMSetLoc(128, 128);
    WMAddFuncButton("leave", &unmapsubmenu, NULL, 0, 0);
    WMDisplay();
    WMAppMainLoop();

    return 0;
}
