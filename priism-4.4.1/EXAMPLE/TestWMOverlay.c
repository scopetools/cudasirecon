/*
 * This program just provides a sampler of how WMSetOverlayUse interacts
 * with other widget creation functions in the WM library.  As a side effect
 * it also gives a brief sample of how to use the various pulldown calls.
 *
 * The overlay planes provide an area to draw without damaging the contents
 * of other planes.  This can be used to avoid expensive redraws.  The problem
 * with the  overlay planes is that they have limited colors (4 - 16 on most
 * Indys and Indigo2s, 256 on O2s and Octanes) and are also used for resizing
 * and moving effects by the window manager.  The implications are that
 * pulldown menus (which are transient, cannot be resized or moved, and have
 * few colors) are generally ok to put in the overlays and this can be a win
 * if the pulldown would normally cause an expensive redraw when it is
 * deactivated by a user.  Putting whole dialogs in the overlay planes
 * does not give a good visual appearance on the systems with limited colors
 * and those dialogs have a different appearance when resized or moved
 * than users would normally expect.
 *
 * By default, nothing is put in the overlay planes.  Call WMSetOverlayUse
 * to change this behavior for all widgets that are created after the call
 * and before the next call to WMSetOverlayUse.  Passing WM_OVERLAY_YES to
 * WMSetOverlayUse specifies that  pulldowns and dialogs should be placed in
 * the overlay planes (whether or not they actually are placed there depends
 * on whether the hardware supports it or not). Passing WM_OVERLAY_NO to
 * WMSetOverlayUse specifies that pulldowns will appear in the same planes as
 * their parents and that dialogs will be placed in the normal planes.  The
 * widget creation calls affected by WMSetOverlayUse are the pulldown calls,
 * WMAddOptionMenu, WMAddPulldown, WMAddPulldownMenu, WMAddNestedPulldown,
 * and WMAddPullRight, and the dialog calls, WMInit, WMInitSubMenu,
 * WMConfirm, WMConfirmError, WMShowInfo, and WMPostInfo (WMGetShellWidget is
 * not affected).
 *
 * Calls to WMOverlayDepth can be used to find the depth of the overlays that
 * the WM library is using.  Since the number of colors available is
 * 2^depth - 1, you can use the return value to determine whether or not to
 * put a dialog or pulldown in the overlay planes at runtime.
 */
    
#include "WMInclude.h"
#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef TRUE
#define TRUE (1 == 1)
#undef FALSE
#define FALSE (1 != 1)


static char moduleName[]  = "TestWMOverlay";

typedef struct DescriptionDialogs {
    int optionSelected;
    Widget dialog;
    Widget dialogInOverlay;
} DescriptionDialogs;


static int HandleFileMenu(int* pOptionSelected);
static int CloseSubDialog(int* pSubdialogID);
static void CreateSubDialog(int putInOverlayPlane);
static int HandleDialogMenu(int* pOptionSelected);
static int HandleHelpMenu(DescriptionDialogs* pDialogInfo);
static void CreateMainDialog(int putInOverlayPlane);


int
main(int argc, char* argv[]) {
    int usageOk;
    int mainDialogInOverlay = FALSE;

    if (argc == 2) {
	if (strcmp("-o", argv[1]) == 0) {
	    usageOk = TRUE;
	    mainDialogInOverlay = TRUE;
	} else {
	    usageOk = FALSE;
	}
    } else {
	usageOk = argc < 2;
    }

    if (usageOk) {
	CreateMainDialog(mainDialogInOverlay);
	WMAppMainLoop();
    } else {
	(void) fprintf(
	    stderr,
	    "usage: %s [-o]\n     -o  put main dialog in overlay plane\n",
	    moduleName
	);
    }

    return (usageOk) ? EXIT_SUCCESS : EXIT_FAILURE;
}


static
void
CreateMainDialog(int putInOverlayPlane)
{
    char* fileOptions[] = { "Exit" };
    static int fileOptionSelected = 0;
    char* dialogOptions[] = {
	"Dialog", "Dialog + Overlay",
	"Confirm Dialog", "Confirm Dialog + Overlay",
	"Error Dialog", "Error Dialog + Overlay",
	"Info Dialog", "Info Dialog + Overlay"
    };
    static int dialogOptionSelected = 0;
    char* helpOptions[] = {
	"Description", "Hide Description",
	"Description + Overlay", "Hide Description + Overlay"
    };
    static DescriptionDialogs helpParams = {0, NULL, NULL};
    Widget pulldownBase;

    WMSetOverlayUse((putInOverlayPlane) ? WM_OVERLAY_YES : WM_OVERLAY_NO);
    WMInit(moduleName);
    WMSetOverlayUse(WM_OVERLAY_NO);
    /*
     * This set of pulldowns always appears in the same planes as
     *its parent.
     */
    pulldownBase = WMAddPulldownMenu(
	NULL,
	"File",
	fileOptions,
	sizeof(fileOptions) / sizeof(fileOptions[0]),
	&fileOptionSelected,
	HandleFileMenu,
	&fileOptionSelected,
	FALSE,
	0
    );
    WMAddPulldownMenu(
	XtParent(pulldownBase),
	"Dialog",
	dialogOptions,
	sizeof(dialogOptions) / sizeof(dialogOptions[0]),
	&dialogOptionSelected,
	HandleDialogMenu,
	&dialogOptionSelected,
	FALSE,
	0
    );
    WMAddPulldownMenu(
	XtParent(pulldownBase),
	"Help",
	helpOptions,
	sizeof(helpOptions) / sizeof(helpOptions[0]),
	&helpParams.optionSelected,
	HandleHelpMenu,
	&helpParams,
	FALSE,
	0
    );
    WMDisplay();

    return;
}


/*
 * Creates a subdialog to demonstrate how the various pulldowns work
 * in combination with calls to WMSetOverlayUse.
 */
static
void
CreateSubDialog(int putInOverlayPlane)
{
    static int subdialogID = 0;
    static int overlaySubdialogID = 0;
    int* pDialogID;

    pDialogID = (putInOverlayPlane) ? &overlaySubdialogID : &subdialogID;

    if (*pDialogID == 0) {
	/*
         * The item at the end is for the menu title when passed to
	 * WMAddPulldown (see the comment with WMSetPulldownShowChoice below.
	 */
	char* menuItems[] = {
	    "Item 1", "Item 2", "Item 3", "Menu"
	};
	char* overlayMenuItems[] = {
	    "Item 1", "Item 2", "Item 3", "OverlayMenu"
	};
	/*
	 * Used in calls to WMAddPullRight so the items in the menu don't lead
	 * to more menus.
	 */
	int atTerminus[] = { 0, 0, 0 };

	char* submenuItems[] = {
	    "Sub", "SubOverlay"
	};
	/*
         * Used in calls to WMAddNestedPulldown so the items in the menu lead
	 * to more menus.
	 */
	int atStart[] = { 1, 1 };

	static int pulldownChoice = 0;
	Widget pulldownBase;
	Widget lastLevelPulldown;

	/*
	 * With this call, the widgets created by WMAddPulldown show a constant
	 * title (specified as the last element in the array of items); using
	 * a non-zero argument would cause the title to display the name of the
	 * item last selected.
	 */
	WMSetPulldownShowChoice(0);
	WMSetOverlayUse((putInOverlayPlane) ? WM_OVERLAY_YES : WM_OVERLAY_NO);
	*pDialogID = WMInitSubMenu("Pulldown Tests");

	/*  This pulldown will appear in the same planes as its parent. */
	WMSetOverlayUse(WM_OVERLAY_NO);
	WMAddPulldown(
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);

	/*  This pulldown always appears in the overlay planes. */
	WMAddPulldown(
	    overlayMenuItems,
	    sizeof(overlayMenuItems) / sizeof(overlayMenuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();

	/*  This pulldown will appear in the same planes as its parent. */
	WMSetOverlayUse(WM_OVERLAY_NO);
	WMAddText("Option", 8, 0); 
	WMAddOptionMenu(
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);

	/* This pulldown always appears in the overlay planes. */
	WMAddText("OverlayOption", 16, 0); 
	WMAddOptionMenu(
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();

	/*
	 * This combination, first pulldown in the menu bar uses normal planes
	 * (assuming the parent dialog is in the normal planes) while the
	 * second is in the overlay planes causes a BadMatch X error at
	 * runtime on SGIs (IRIX Motif, as an optimization, uses the
         * same shell widget for the pulldowns in a menu bar).
	 */
	/*
	WMSetOverlayUse(WM_OVERLAY_NO);
	pulldownBase = WMAddPulldownMenu(
	    NULL,
	    "Menu",
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPulldownMenu(
	    XtParent(pulldownBase),
	    "OverlayMenu",
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();
	*/

	/*
	 * Only the first WMSetOverlayUse has any effect.  All pulldowns
	 * from the menubar end up in the overlay planes.
	 */
	WMSetOverlayUse(WM_OVERLAY_YES);
	pulldownBase = WMAddPulldownMenu(
	    NULL,
	    "OverlayMenu",
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();
	WMSetOverlayUse(WM_OVERLAY_NO);
	WMAddPulldownMenu(
	    XtParent(pulldownBase),
	    "Menu",
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();

	/*
         * When the subdialog is not in the overlay planes, the first pullright
	 * is in the normal planes and the second is in the overlay planes.
         * When the subdialog is in the overlay planes, both pullrights are in
         * the overlay planes.  This is not of much practical use of course.
         * Typical usage would be to just have WMSetOverlayUse(WM_OVERLAY_YES)
         * before the parent WMAddNestedPulldown.
	 *
	 * This combination causes bad match errors with OpenMotif 2.1.30.
	 */
	/*
	WMSetOverlayUse(WM_OVERLAY_NO);
	pulldownBase = WMAddNestedPulldown(
	    NULL,
	    "Menu",
	    atStart,
	    submenuItems,
	    sizeof(submenuItems) / sizeof(submenuItems[0]),
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAddPullRight(
	    pulldownBase,
	    0,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPullRight(
	    pulldownBase,
	    1,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();
	*/

	/*
	 * Both pullrights are always in the overlay planes.  The extra calls
	 * to WMSetOverlayUse after the first have no effect on the appearance.
	 */
	WMSetOverlayUse(WM_OVERLAY_YES);
	pulldownBase = WMAddNestedPulldown(
	    NULL,
	    "OverlayMenu",
	    atStart,
	    submenuItems,
	    sizeof(submenuItems) / sizeof(submenuItems[0]),
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_NO);
	WMAddPullRight(
	    pulldownBase,
	    0,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPullRight(
	    pulldownBase,
	    1,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();

	/*
	 * Causes BadMatch errors at runtime on SGIs.
	 */
	/*
	WMSetOverlayUse(WM_OVERLAY_NO);
	pulldownBase = WMAddNestedPulldown(
	    NULL,
	    "Menu",
	    atStart,
	    submenuItems,
	    sizeof(submenuItems) / sizeof(submenuItems[0]),
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAddPullRight(
	    pulldownBase,
	    0,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPullRight(
	    pulldownBase,
	    1,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	lastLevelPulldown = WMAddNestedPulldown(
	    XtParent(pulldownBase),
	    "OverlayMenu",
	    atStart,
	    submenuItems,
	    sizeof(submenuItems) / sizeof(submenuItems[0]),
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_NO);
	WMAddPullRight(
	    lastLevelPulldown,
	    0,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems)/sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPullRight(
	    lastLevelPulldown,
	    1,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();
	*/

	/*
	 * All the pulldowns and pullrights in this hierarchy end up in
	 * the overlay planes.  Only the first call to WMSetOverlayUse
	 * has an effect on the appearance.
	 */
	WMSetOverlayUse(WM_OVERLAY_YES);
	pulldownBase = WMAddNestedPulldown(
	    NULL,
	    "OverlayMenu",
	    atStart,
	    submenuItems,
	    sizeof(submenuItems) / sizeof(submenuItems[0]),
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_NO);
	WMAddPullRight(
	    pulldownBase,
	    0,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPullRight(
	    pulldownBase,
	    1,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_NO);
	lastLevelPulldown = WMAddNestedPulldown(
	    XtParent(pulldownBase),
	    "Menu",
	    atStart,
	    submenuItems,
	    sizeof(submenuItems) / sizeof(submenuItems[0]),
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAddPullRight(
	    lastLevelPulldown,
	    0,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMSetOverlayUse(WM_OVERLAY_YES);
	WMAddPullRight(
	    lastLevelPulldown,
	    1,
	    atTerminus,
	    menuItems,
	    sizeof(menuItems) / sizeof(menuItems[0]) - 1,
	    &pulldownChoice,
	    NULL,
	    NULL,
	    0,
	    0
	);
	WMAttachRightSide();
	WMNewRow();

	/*
	 * The function button doesn't have it's own shell so its always in the
	 * sames planes as its parent irregardless of prior calls to
         * WMSetOverlayUse.
	 */
	WMAddFuncButton("Close", CloseSubDialog, pDialogID, 0, 0);
	WMAttachRightSide();
    }

    WMDisplaySubMenu(*pDialogID);

    return;
}


static
int
CloseSubDialog(int* pSubdialogID)
{
    WMUndisplaySubMenu(*pSubdialogID);
    return 0;
}


/*
 * Demonstrate creating each type of dialog either in the normal planes
 * or overlay planes.
 */
static
int
HandleDialogMenu(int* pOptionSelected)
{
    if (pOptionSelected == NULL) {
	(void) fprintf(
	    stderr, "%s: Null pointer in HandleDialogMenu\n", moduleName
	);
    } else {
	int inOverlay = TRUE;

	switch (*pOptionSelected) {

	case 0:
	    inOverlay = FALSE;
	    /* fall through */
	case 1:
	    CreateSubDialog(inOverlay);
	    break;

	case 2:
	    inOverlay = FALSE;
	    /* fall through */
	case 3:
	    WMSetOverlayUse(inOverlay);
	    WMConfirm("How's this?");
	    break;

	case 4:
	    inOverlay = FALSE;
	    /* fall through */
	case 5:
	    WMSetOverlayUse(inOverlay);
	    WMConfirmError("And how's this?");
	    break;

	case 6:
	    inOverlay = FALSE;
	    /* fall through */
	case 7:
	    WMSetOverlayUse(inOverlay);
	    /*  Keep it in place for at most 5 seconds. */
	    WMPostInfo("And what about this?", 5);
	    break;

	default:
	    (void) fprintf(
		stderr,
		"%s: Unrecognized option in HandleDialogMenu\n",
		moduleName
	    );
	    break;
	}
    }
    return 0;
}


static
int
HandleFileMenu(int* pOptionSelected)
{
    if (pOptionSelected == NULL) {
	(void) fprintf(
	    stderr, "%s: Null pointer in HandleFileMenu\n", moduleName
	);
    } else {
	switch (*pOptionSelected) {

	case 0:
	    exit(EXIT_SUCCESS);
	    break;

	default:
	    (void) fprintf(
		stderr, "%s: Unrecognized option in HandleFileMenu\n",
		moduleName
	    );
	    break;
	}
    }
    return 0;
}


static
int
HandleHelpMenu(DescriptionDialogs* pDialogInfo)
{
    if (pDialogInfo == NULL) {
	(void) fprintf(
	    stderr, "%s: Null pointer in HandleHelpMenu\n", moduleName
	);
    } else {
	char message[] = "Test of WM library overlay functions";

	switch (pDialogInfo->optionSelected) {
	case 0:
	    if (pDialogInfo->dialog != NULL) {
		WMRemoveInfo(pDialogInfo->dialog);
	    }
	    WMSetOverlayUse(WM_OVERLAY_NO);
	    pDialogInfo->dialog = WMShowInfo(message);
	    break;

	case 1:
	    if (pDialogInfo->dialog != NULL) {
		WMRemoveInfo(pDialogInfo->dialog);
		pDialogInfo->dialog = NULL;
	    }
	    break;

	case 2:
	    if (pDialogInfo->dialogInOverlay != NULL) {
		WMRemoveInfo(pDialogInfo->dialogInOverlay);
	    }
	    WMSetOverlayUse(WM_OVERLAY_YES);
	    pDialogInfo->dialogInOverlay = WMShowInfo(message);
	    break;

	case 3:
	    if (pDialogInfo->dialogInOverlay != NULL) {
		WMRemoveInfo(pDialogInfo->dialogInOverlay);
		pDialogInfo->dialogInOverlay = NULL;
	    }
	    break;

	default:
	    (void) fprintf(
		stderr,
		"%s: Unrecognized option in HandleHelpMenu\n",
		moduleName
	    );
	    break;
	}
    }
    return 0;
}
