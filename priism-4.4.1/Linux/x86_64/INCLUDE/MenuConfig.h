#ifndef UCSF_MSG_MENUCONFIG_H
#define UCSF_MSG_MENUCONFIG_H

/*
 * Routines for building a menu bar whose entries consist of a mixture
 * of submenus, actions for launching external applications, application-
 * defined actions, and help items.  A configuration file can be read to
 * define some or all of the submenus, external application actions, and
 * help items.
 *
 * The procedure for using these routines is:
 *
 * 1) Initialize the menu description with InitializeMenuDescription.
 *
 * 2) (optional) Add application-defined items using AppendMenuEntries
 *    and FindSubmenu (to locate the submenu to which the entries will be
 *    added).
 *
 * 3) (optional) Read from a configuration file and add the entries it
 *    describes using AppendMenusFromFile.  When entries are added to a
 *    menu that already has entries from step 2, the new entries appear
 *    after the ones already present.
 *
 * 4) (optional) Add more application-defined items using the same process
 *    as step 2.  When adding to menus which already contain items, the
 *    new items appear after the existing ones.
 *
 * 5) (optional) Repeat steps 3 and 4 as desired.
 *
 * 6) (optional) Using the entries present in the menu description so far,
 *    build a similar hierarchy of items in the Help menu (if the Help menu
 *    exists, this hierarchy of items appears after the existing items.
 *    GenerateHelpItems is the call used to do this.
 *
 * 7) (optional) Repeat steps 3 and 4 to add items that you didn't want in
 *    the menus when AppendHelpMenus did its work (best used for adding
 *    items to the help menu after anything added in step 6).
 *
 * 8) Compile the menu description and create the widgets (with
 *    a series of WMAddNestedPulldown and WMAddPullright calls) for the menus
 *    using CompileMenu.  Prior to doing this you must have initialialized
 *    the WM library.  
 *
 * Any time after step 2 and before step 8, you can mark any submenu that
 * you'll want to modify after compilation using MarkAlterableSubmenu.  Then
 * after step 8 you can modify the labels (using AlterSubmenuLabel) or
 * sensitivities of the elements of the submenu and then use UpdateSubmenu
 * to commit the changes.
 *
 * At any point after step 8, you can free the original menu description
 * using FreeMenuDescription.  The compiled menu information can only
 * be freed safely using FreeCompiledMenu after the dialog in which the menu
 * is embedded is destroyed.
 */


#include <X11/Intrinsic.h> /* for Widget */
#include <stddef.h> /* for size_t */


struct MenuDescription;
typedef struct MenuDescription MenuDescription;

struct Submenu;
typedef struct Submenu Submenu;

struct CompiledMenu;
typedef struct CompiledMenu CompiledMenu;


/*
 * All allocations are done using the function referenced by allocator.
 * deallocator is used to reverse allocator.  For either of these function
 * pointers, you can specify NULL to get a default allocation/deallocation
 * routine.
 *
 * NULL is returned if there was a failure in the initial allocation of
 * space for the menu description.  A successful return value is a non-NULL
 * pointer.
 */
MenuDescription* InitializeMenuDescription(
    void* (*allocator)(size_t byteCount),
    void (*deallocator)(void* ptr)
);

/*
 * No effect on a NULL pointer.
 */
void FreeMenuDescription(MenuDescription* pDescr);


/*
 * Searches the contents of a submenu for another submenu with the label
 * labelToFind.  If no submenu is specified (pSearchStart is NULL),
 * the labels of the top level menus are searched.
 *
 * NULL is returned if the menu was not found or pDescr or pSearchStart
 * could be identified as invalid pointers.
 */
Submenu* FindSubmenu(
    MenuDescription* pDescr, Submenu* pSearchStart, const char* labelToFind
);


/*
 * If the first character in filename is not a '/', then filename is used
 * as the last part of the path and $(HOME)/.iveprefs is tried as the first
 * part of the path and, if that does not locate a file, then the configuration
 * directory in the IVE distribution is used as the first part of the path.
 * Otherwise the filename is used without modification.
 *
 * 0 is returned upon successful completion, otherwise a non-zero value
 * is returned and the description is left in a compilable but likely
 * incomplete state.
 *
 * Error messages from parsing the file are written to standard error.
 * Standard input is reopened to parse the input file.
 */
int AppendMenusFromFile(MenuDescription* pDescr, const char* filename);

/*
 * Appends the first count items from the pEntries array to the menu
 * referenced by pSubmenu or to the top-level menu if pSubmenu is NULL.
 *
 * 0 is returned upon successful completion, otherwise a non-zero value
 * is returned and the description is left in a compilable but likely
 * incomplete state.
 */
typedef enum MenuEntryOption {
    MenuEntryExecItem = 0,
    MenuEntryUserItem = 1,
    MenuEntryHelpItem = 2,
    MenuEntrySeparator = 3,
    MenuEntrySubmenu = 4
} MenuEntryOption;

typedef struct ExecItemSpec {
    char* label;
    char* pathToExecutable;
    char* helpLookup;
} ExecItemSpec;

typedef struct UserItemSpec {
    char* label;
    char* helpLookup;
    void (*callback)(void*);
    void* param;
} UserItemSpec;

typedef struct HelpItemSpec {
    char* label;
    char* helpLookup;
} HelpItemSpec;

typedef struct SubmenuSpec {
    char* label;
} SubmenuSpec;

typedef struct MenuEntry {
    union {
	ExecItemSpec execItem;
	UserItemSpec userItem;
	HelpItemSpec helpItem;
	SubmenuSpec submenu;
    } specif;
    MenuEntryOption type;
} MenuEntry;

int AppendMenuEntries
(
    MenuDescription* pDescr,
    Submenu* pSubmenu,
    const MenuEntry* pEntries,
    int count
);

/*
 * Using the current contents of the menus (excluding the Help menu), appends
 * to the help menu (creating it if necessary) a hierarchy of links to
 * on-line help mirroring the menu hierarchy.
 *
 * 0 is returned upon successful completion, otherwise a non-zero value
 * is returned and the description is left in a compilable but likely
 * incomplete state.
*/
int GenerateHelpItems(MenuDescription* pDescr);


/*
 * Returns a NULL pointer if a memory allocation error occurred during the
 * compile or pDescr could be identified as an invalid pointer.
 * Otherwise a non-NULL pointer is returned.
 *
 * If firstMenuRtn is non-NULL and the compilation is successful, the
 * location it references is set to the id of the first top-level menu in the
 * menu system.
 */
CompiledMenu* CompileMenu
(
 Widget* firstMenuRtn,
 const MenuDescription* pDescr,
 void* (*allocator)(size_t byteCount),
 void (*deallocator)(void* ptr),
 void (*alternateExecCallback)(void* pathAsVoidPtr)
);

/*
 * No effect on a NULL ptr.
 */
void FreeCompiledMenu(CompiledMenu* ptr);


typedef struct SubmenuAlterables {
    /* Not alterable but for calls to WMChangePulldown. */
    Widget submenu;
    int* sensitivities;
    /* Use AlterSubmenuLabel to change an element of this array. */
    char** labels;
    /* Not alterable, used by AlterSubmenuLabel. */
    CompiledMenu* pCompiledMenu;
    /* Not alterable but for calls to WMChangePulldown. */
    int count;
} SubmenuAlterables;

/*
 * Marks a submenu as alterable (or not alterable) after the compilation step.
 *
 * To mark a submenu as alterable, pass a non-NULL address for the first
 * parameter and the address of the submenu as the second parameter.
 * When the menu is successfully compiled, the value pointed to by the
 * first parameter will be set to the location of a structure storing the
 * alterable parameters for the menu.  After altering the parameters, use
 * them in a call to WMChangePulldown to change the physical appearance of
 * the pulldown.  When FreeCompiledMenu is called, the memory for this
 * structure is freed, so accesses subsequent to a FreeCompiledMenu call
 * would give undefined behavior.
 *
 * To mark the menu as non-alterable (the default), pass NULL as the first
 * parameter and the address of the submenu as the second parameter.
 */
void MarkAlterableSubmenu(SubmenuAlterables** ppMemento, Submenu* pSubmenu);

/*
 * To safely change the ith label in an alterable submenu use this function.
 * 0 is returned upon successful completion; a non-zero value is returned
 * upon failure and the submenu alterable data is left unaltered.
 */
int AlterSubmenuLabel(SubmenuAlterables* ptr, int index, const char* newLabel);


/*
 * Since WMDeleteWidget, WMPasteWidget, WMUnpasteWidget only operate
 * on a single submenu, it is convenient to know all the menus in a menu
 * bar.
 *
 * The call below sets the location where the location of a NULL terminated
 * list of the toplevel menu's widget IDs should be stored.  This call
 * must be made before CompileMenu which generates the list.  The list is
 * freed when FreeCompiledMenu is called.
 *
 * 0 is returned on success and 1 is returned on failure (ppList was
 * not NULL or a valid pointer or pDescr was NULL or an invalid pointer).
 */
int SetToplevelListLoc(MenuDescription* pDescr, Widget** ppList);

#endif /* include guard */
