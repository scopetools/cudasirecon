#ifndef UCSF_MSG_WMFUNC_H
#define UCSF_MSG_WMFUNC_H


#include <GL/glx.h> /* for GLXContext */
#include <X11/Intrinsic.h> /* Widget, XtAppContext, XtPointer,
			      XtTimerCallbackProc */


/*
 * Traditionally, the C prototypes for the WM functions have not specified
 * what the arguments are for the callback functions, and in C++
 * they have either used int*, char*, or void in the prototype.  For the
 * record, here are the definitions that correspond to how the functions are
 * invoked by the library (userParam is the int* or int (depending on which
 * callback is involved) that is passed when the callback was registered).
 * The only callback return value that is used by the library is that for the
 * WMDoButtonCB (if < 0 nothing is appended to the script and nothing is
 * submitted).
 *
 * typedef int (*WMCallback)(int* userParam);
 * typedef int (*WMCharFieldCB)(int* userParam);
 * typedef int (*WMDoButtonCB)(int* userParam, char* userCommandString);
 * typedef int (*WMScrolledChoicesCB)(int selection, int* userParam);
 * typedef void (*WMGLCallback)(Widget drawarea, XtPointer userParam, XtPointer pGLDrawAreaCallbackStruct);
 * typedef void (*WMExitCB)(void);
 * typedef int (*WMSubExitCB)(void)
 * typedef int (*WMDisplayChangeCB)(int userParam, XEvent* event);
 * typedef int (*WMDisplayChangeFortranCB)(int* userParam, int* eventType);
 * typedef int (*WMEventCB)(int userParam, XEvent* event);
 * typedef int (*WMClientMessageCB)(Widget form, int userParam, XEvent* event);
 *
 * The definitions below are how these functions have been traditionally
 * prototyped (with the exception of the C++ prototypes for WMAddOptionMenu,
 * WMAddPulldownMenu, WMAddPullRight, WMAddNestedPulldown which now match
 * how they are used by the library and declared in C (with the exception of
 * the first, these functions are newer and none were used in C++ code in
 * Priism as of December 1999 so the backward compatibility problem seemed
 * less important than consistency)).
 */
#ifdef __cplusplus
extern "C"
{
typedef int (*WMCallback)(int*);
typedef int (*WMCharFieldCB)(char*);
typedef int (*WMDoButtonCB)(void);
typedef int (*WMScrolledChoicesCB)(int*);
typedef int (*WMGLCallback)(int*);
typedef int (*WMExitCB)(int*);
typedef int (*WMSubExitCB)(int*);
typedef int (*WMDisplayChangeCB)(void);
typedef int (*WMDisplayChangeFortranCB)(void);
typedef int (*WMEventCB)(void);
typedef int (*WMClientMessageCB)(void);
#else
typedef int (*WMCallback)();
typedef int (*WMCharFieldCB)();
typedef int (*WMDoButtonCB)();
typedef int (*WMScrolledChoicesCB)();
typedef int (*WMGLCallback)();
typedef int (*WMExitCB)();
typedef int (*WMSubExitCB)();
typedef int (*WMDisplayChangeCB)();
typedef int (*WMDisplayChangeFortranCB)();
typedef int (*WMEventCB)();
typedef int (*WMClientMessageCB)();
#endif


XtAppContext WMRetAppContext(void);
void WMInitApp(void);
Widget WMGetShellWidget(const char* title);
void WMSync(void);
Widget WMInit(const char* title);
void WMForceUpdate(Widget w);
void WMProcessEvent(void);
int WMInitSubMenu(const char* title);
int WMDisplaySubMenu(int w);
int WMUndisplaySubMenu(int w);
XEvent* WMRetLastButEvent(void);
int WMGetLoc(int* pX, int* pY, int menuNumber);
int WMGetMenuSize(int* width, int* height, int menu_number);
int WMGetScreenSize(int* width, int* height);
int WMSetLoc(int x, int y);
int WMDeleteWidget(Widget w,int resize_flag);
int WMInsertWidget(Widget w,int method);
int WMEndInsert(int resize_flag);
int WMUnpasteWidget(Widget w);
int WMPasteWidget(Widget w);
int WMEnableField(Widget w);
int WMDisableField(Widget w);
void WMEnableGroup(int group);
void WMDisableGroup(int group);
int WMSensitiveGroup(int group);
int WMInsensitiveGroup(int group);
int WMUpdateGroup(int flag);
int WMUpdateField(Widget w);
int WMUpdateFuncList(int ilist, const char* const* funcname,
		     const char* const* helpindex, int nfunc, int update);
Widget WMAddArrowButton(int dir, WMCallback callback, int* param);
void WMInitChainCMD(void);
void WMStopChainCMD(void);
int WMSetOverlayUse(int overlayState);
int WMOverlayDepth(void);
Widget WMAddDoButton(WMDoButtonCB build_com, int* param, const char* qcom,
		     char* command);
Widget WMAddSaveButton(const char* filename);
Widget WMNewRow(void);
Widget WMAddSeparator(void);
Widget WMAddScrolledChoices(char* strings[], int* num_items_ptr,
			    int num_visible, WMScrolledChoicesCB callback,
			    int* params, int group);
int WMChangeScrolledChoicesSelection(Widget w, int new_selection);
Widget WMAddScrolledText(
             const char* text, const int* nitem, int nvis, int group);
Widget WMAddFuncList(const char* const* funcname, const char* const* helpindex,
             int nfunc, int width, int height, int* ifunc, WMCallback callback,
	     int* param, int update);
void WMSetGetFileSrc(const char* dirkey, const char* patkey);
int WMGetFile(char* filename, char* dir, char* pattern, char* filter);
Widget WMAddGetFile(const char* label, char* dir, char* pattern, int maxlen,
                    int dislen, char* filename, WMCallback callback,
		    int* param, int update, int group);
Widget WMAddText(char* text, int maxlen, int group);
Widget WMAddGLwMDrawWidget(int width, int height,
           WMGLCallback expose_callback, int* ex_param,
           WMGLCallback input_callback, int* input_param,
           GLXContext* context);
Widget WMAddDrawingArea(int width, int height);
int WMChangeText(Widget widget, const char* new_string);
int WMChangePulldown(Widget widget, const char* const* new_labels,
		      const int sens_list[], int nlabel);
Widget WMAddFloatField(float* fptr, int num, int dislen, int maxlen,
                       float* min, float* max, int* decimal,
		        WMCallback callback, int* param, int update, int group);
Widget WMAddIntField(int* iptr, int num, int dislen, int maxlen,
                     int* min, int* max, int* factor,
		      WMCallback callback, int* param, int update, int group);
Widget WMAddCharField(char* text, int maxlen, int dislen,
		       WMCharFieldCB callback, int* param,
		       int update, int group);
int WMSetPulldownShowChoice(int i);
int WMGetPulldownShowChoice(void);
Widget WMAddPulldown(const char* const* label, int nlabel, int* ichoice,
                     WMCallback callback, int* param, int update, int group);
Widget WMAddPulldownMenu(Widget ref, const char* menu_name,
			  const char* const* label,int nlabel, int* ichoice,
			  WMCallback callback, XtPointer param,
                         int update, int group);
Widget WMAddPullRight(Widget ref, int element, const int add_level[],
		       const char* const* label, int nlabel, int* ichoice,
		       WMCallback callback, XtPointer param,
		       int update, int group);
Widget WMAddNestedPulldown(Widget ref, const char* menu_name,
			    const int addlevel[], const char* const* label,
			    int nlabel, int* ichoice, WMCallback callback,
			    int* param, int update, int group);
Widget WMAddOptionMenu(const char* const* label, int nlabel, int* ichoice,
                       WMCallback callback, XtPointer param, int update,
		        int group);
int WMSliderShowValue(int flag);
void WMAttachRightSide(void);
Widget WMAddHorSlider(const char* label, int width, int min, int max,
		       int* ivar, WMCallback scallback, int *param,
		       WMCallback dcallback, int *dparam, int update,
		       int group);
Widget WMAddToggleButton(const char* label, int* ivar, WMCallback callback,
                         int* param, int update, int group);
Widget WMAddOnOffStatus(const char* label, int* ivar, int group);
Widget WMAddStatusBar(const char* label, int width, int maxval, int* ivar,
		       int update, int group);
Widget WMAddFuncButton(char* label, WMCallback callback, int* param,
                       int update, int group);
Widget WMAddInfoButton(const char* label, const char* keyword);
int WMSetOffset(int ival1, int ival2, int ival3, int ival4);
int WMDisplay(void);
void WMEnableIWLEvent(void);
void WMDisableIWLEvent(void);
int WMProcEvent(int wait);
int WMAppMainLoop(void);
int WMConfirm(const char* message);
int WMConfirmError(const char* message);
int WMPostInfo(const char* message, int second);
int WMNotExecutable(const char* name);
Widget WMShowInfo(const char* message);
int WMRemoveInfo(Widget wigt);
int WMRegMouseButton(int istream, int button_mask);
int WMUnRegMouseButton(int istream, int button_mask);
int WMAddEventHandler(int istream, long mask, WMEventCB function, int argu);
int WMCancelEventHandler(int istream, long mask);
int WMProcExtEvent(int event_type, int (*function)(), int argu);
int WMCancelExtEvent(int event_type);
int WMProcDisplayChange(int istream, WMDisplayChangeCB function, int argu);
int WMCancelDisplayChange(int istream);
int WMProcClientMessage(WMClientMessageCB function, int argu);
int WMCancelClientMessage(void);
void WMDisplayHelp(const char* keyword);
int WMRingBell(int ivol, float len);
int WMGetLastButNum(void);
int WMNoTitleBar(void);
int WMOverrideRedirect( int flag );
int WMOpenSclFile(void);
int WMGetSclRec(char*, SCALE*);
int WMSaveSclRec(char*, SCALE*);
void WMCloseSclFile(void);
void WMSetExitFunction(WMExitCB callback);
void WMSetSubExitFunc(int sm, WMSubExitCB callback);
int WMChangeTitle(Widget widget, const char* new_string);
int WMStartRemoteProgram(const char* machine_name, const char* full_command);
int WMAddInputEventHandler(int (*function)(), int ref, int argu);
int WMCancelInputEventHandler(int handler_ref);
int WMRemoteSendData(const void* data, int handler_code, int size);
int WMRemoteGetData(void* data, int size);
int WMRemoteGetBytes(void* data, int size);
int WMRemoteGetInt2(void* data, int nget);
int WMRemoteGetInt4(void* data, int nget);
int WMOpenNoteFile(const char* filename);
int WMCloseNoteFile(const char* filename);
int WMUpdateNoteFile(const char* filename);
int WMSetIntNote(const char* label, int nvalue, const int array[]);
int WMSetFloatNote(const char* label, int nvalue, const float array[]);
int WMSetStringNote(const char* label, int nvalue, const char* const* str);
int WMGetIntNote(const char* label, int* nvalue, int array[]);
int WMGetFloatNote(const char* label, int* nvalue, float array[]);
int WMGetStringNote(const char* label, int* nvalue, char** str);
int WMGetStringNoteLength(const char* label, int* nvalue, int* length);

/*
 * Only used internally?
 */
Window GetXWindowID(void);


#ifdef __cplusplus
} /* close extern "C" */
#endif

#endif /* include guard */
