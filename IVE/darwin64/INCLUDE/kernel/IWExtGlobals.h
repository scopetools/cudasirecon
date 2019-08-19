#ifndef UCSF_MSG_IWEXTGLOBALS_H
#define UCSF_MSG_IWEXTGLOBALS_H

/*-----------------------------------------------------------------------------
    Copyright (C) 1993
    Macromolecular Structure Group of Biochemistry Dept. at University of
    California at San Francisco.
    These coded instructions, statements, and computer programs comprise
    unpublished proprietary information of the Macromolecular Structure Group of
    the Biochemistry Department at University of California at San Francisco, 
    and are protected by Federal copyright law.  They may not be disclosed
    to third parties, copied or duplicated in any form - in whole or in part -
    without the prior written consent of Macromolecular Structure Group of
    Biochemistry Department at University of California at San Francisco.
-------------------------------------------------------------------------*/

/*
 * Files including this file should include ive_standards.h before any headers
 * (see note in kernel/IWConstants.h).
 */
#include "kernel/IWConstants.h"  /* IW_SHARED_MEM_STRUCT */
#include "ive_shm.h"      /* IVEArena */
#include <X11/Xlib.h>     /* Display */

extern IW_SHARED_MEM_STRUCT* IW_iwl_ptr;
extern IVEArena IW_arena;
extern char* IW_color_names[];
extern Display* IW_appl_display;

#endif /* include guard */
