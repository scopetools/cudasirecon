#ifndef UCSF_MSG_WMINCLUDE_H
#define UCSF_MSG_WMINCLUDE_H

/*-----------------------------------------------------------------------------
    Copyright (C) 1993
    Macromolecular Structure Group of Biochemistry Dept. at University of
    California at San Francisco
    These coded instructions, statements, and computer programs comprise
    unpublished proprietary information of the Macromolecular Structure Group of
    the Biochemistry Department at University of California at San Francisco,
    and are protected by Federal copyright law.  They may not be disclosed
    to third parties, copied or duplicated in any form - in whole or in part -
    without the prior written consent of Macromolecular Structure Group of
    Biochemistry Department at University of California at San Francisco.
-------------------------------------------------------------------------*/


#include <X11/keysym.h>


typedef struct scale_ {
   char   filename[200];
   int    iwave[5];
   float  scale[5][4];
   int    xst,width,yst,height;
   int    zst, zend, zinc;
   int    tst, tend, tinc;
} SCALE, SCALE_PTR;

#include "WMFunc.h"

#define WM_UP_ARROW              0
#define WM_DOWN_ARROW            1
#define WM_LEFT_ARROW            2
#define WM_RIGHT_ARROW           3

#define WM_OVERLAY_NO            0
#define WM_OVERLAY_YES           1


#endif /* include guard */
