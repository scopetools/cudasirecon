#ifndef UCSF_MSG_IP_H
#define UCSF_MSG_IP_H

#include "ipconstants.h"
#include "IWInclude.h"     /* IW_MAX_WAVE */
#include <X11/Intrinsic.h> /* for Widget */

#ifdef __cplusplus
typedef void (*IPCallback)(int);
extern "C" {

#else
/*
 * For backward compatibility, declare without specifying the argument
 * (internally an int is passed to the callback except for PROC_FUNC functions
 * in 2D processing (they are passed a three element integer array, a three
 * element floating-point array, and an integer) and ENV_CMD functions
 * (which are passed an integer and an array of strings).
 */
typedef void (*IPCallback)();

#endif /* else clause of ifdef __cplusplus */

extern void   IPAddInput(
    int in_stream,
    const char* label,
    int required,
    int win_only,
    IPCallback file_open_hdl,
    int argu
);
extern int    IPAddOutput(
    int out_stream,
    const char* label,
    int required,
    int win_only,
    IPCallback file_open_hdl,
    int argu
);
extern void   IPAllowMode(int out_stream, int mode, int yesno);
extern void   IPAllowModeChange(int out_stream, int yesno);
extern void   IPAppendCommand(const char* text);
extern int    IPCheckInterrupt(void);
extern void   IPDesist(void);
extern void   IPDisplayMessage(const char* msg, int priority);
extern void   IPEnableDoSec(void);
extern void   IPEnableRegionSelection(int in_stream);
extern void   IPEnableResolutionSelection(void);
extern void   IPGetAppendOptions(
    int out_stream, int* p_in_z, int* p_in_w, int* p_in_t, int* p_replace
);
extern char*  IPGetArg(int iarg);
extern int    IPGetArgCount(void);
extern int    IPGetCurRes(int im_stream);
extern void   IPGetCurZWT(int izwt[6]);
extern void   IPGetCurZWTExt(int im_stream, int izwt[3]);
extern float* IPGetDataPtr(int im_stream);
extern int    IPGetDefOutPutMode(void);
extern int    IPGetDefOutPutModeExt(int out_stream);
extern void   IPGetFilename(int istream, char* file);
extern void   IPGetRegionInfo(int nxyzt[4][3], int iwave[IW_MAX_WAVE]);
extern void   IPGetRegionInfoExt(
    int in_stream, int nxyzt[4][3], int iwave[IW_MAX_WAVE]
);
extern void   IPGetRegionWidgets(Widget widgets[IW_MAX_WAVE + 3]);

extern void   IPGetWaveOrder(int in_stream, int* p_nout, int waves[]);
extern int    IPHasGUI(void);
extern int    IPIsAppending(int out_stream);
extern int    IPIsInteractive(void);
extern int    IPIsOpen(int im_stream);
extern int    IPParseIntArg(
    const char* arg, const char* prefix, int max_count, int values[]
);
extern int    IPParseRealArg(
    const char* arg, const char* prefix, int max_count, float values[]
);
extern void   IPPostInfo(const char* string);
extern int    IPPromptInt(
    const char* name, int set_count, int min_count, int max_count, int values[]
);
extern int    IPPromptReal(
    const char* name,
    int set_count,
    int min_count,
    int max_count,
    float values[]
);
extern int    IPPromptString(const char* prompt, char* value, int max_length);
extern void   IPReportSectionStats(
    int out_stream, float smin, float smax, float ssum
);
extern void   IPSetAppendOptions(
    int out_stream, int in_z, int in_w, int in_t, int replace
);
extern void   IPSetAutoScale(int onoff);
extern void   IPSetComputeSize(
    int out_stream, void (*func)(int istrm, int arg, int sizes[5]), int argu
);
extern void   IPSetCusProc(int func_op, IPCallback func, int argu);
extern void   IPSetCustRoutine(int func_op, IPCallback func, int argu);
extern void   IPSetDataPtr(int im_stream, float* array);
extern void   IPSetDefOutPutMode(int mode);
extern void   IPSetDefOutPutModeExt(int out_stream, int mode);
extern void   IPSetDisChgHandler(int im_stream, IPCallback func, int argu);
extern void   IPSetFirstWaveMMM(float fmin, float fmax, float fmean);
extern void   IPSetFirstWaveMMMExt(
    int out_stream, float fmin, float fmax, float fmean
);
extern void   IPSetIncoreProc(void);
extern void   IPSetIncoreProcExt(int out_stream, int partner_stream);
extern void   IPSetMenuLoc(int xpos, int ypos);
extern void   IPSetMenuTitle(const char* title);
extern void   IPSetNoAppend(void);
extern void   IPSetNoAutoIO(void);
extern void   IPSetProcType(int itype);
extern void   IPSetSclParam(int ichan, const float scale[4]);
extern void   IPSetSclParamExt(
    int out_stream, int ichan, const float scale[4]
);

#ifdef __cplusplus
} /* close extern "C" */
#endif /* ifdef __cplusplus */

#endif /* include guard */
