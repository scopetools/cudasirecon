#ifndef UCSF_MSG_IK_H
#define UCSF_MSG_IK_H

/*
 * Copyright(C) 2001
 * Macromolecular Structure Group of the Biochemistry Department at the
 * University of California, San Francisco
 */

#include <sys/types.h>  /* off_t */


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int IWResolution_nres_sec(int wid);
int IWResolution_sec(int wid, int sec, int cur_res);
int IKStartMonitor(int wid);
int IKGetSHMSize(void);

int IKRtDisArea(int wid, int* ix, int* iy, int* width, int* height);
int IKRtWinPos(int wid, int* xst, int* yst, int* width, int* height);
int IKAlWinPos(int wid, int xst, int yst, int width, int height);
int IKRtMulDis(int wid, int* mdx, int* mdy);
int IKAlMulDis(int wid, int mdx, int mdy);
int IKRtDisFmt(int wid);
int IKAlDisFmt(int wid, int iorder);
int IKRtColorMode(int wid);
int IKAlColorMode(int wid, int mode);
int IKRtImgColor(int wid, int colors[]);
int IKAlImgColor(int wid, const int colors[]);
int IKAlGrDisRange(int wid, const int* range);
int IKRtGrDisRange(int wid, int* range);
int IKRtPsdGrColor(int wid);
int IKAlPsdGrColor(int wid, int value);
int IKAlDisWinBdr(int wid, int border_settings);
int IKRtDisWinBdr(int wid, int* p_border_settings);
int IKAlDisWinTool(int wid, int tool_settings);
int IKRtDisWinTool(int wid, int* p_tool_settings);
int IKRtTopImgLoc(int wid, int wave);
int IKRtScratchSize(int wid);
int IKAlScratchSize(int wid, int size);
int IKClrScr(int wid);
int IKAlClrBkg(int wid, int flag);
int IKIsWinMap(int wid);
int IKUnmapWin(int wid);
int IKMapWin(int wid);
int IKRtDisDelay(int wid);
int IKAlDisDelay(int wid, int delay);
int IKAlInterp(int wid, int flag);
int IKRtInterp(int wid);
int IKAlDisImg(int wid, int flag);
int IKRtDisImg(int wid);
int IKAlBufImg(int wid, int flag);
int IKRtBufImg(int wid);
int IKAlNoCopy(int wid, int flag);
int IKRtNoCopy(int wid);
int IKAlFileFormat(int wid, int flag);
int IKSetMemSize(int wid);
int IKAttachFile(int wid, const char* filename);
float IKRtZoom(int wid);
int IKAlZoom(int wid, float zoom);
int IKIsWaveMap(int wid, int wave);
int IKMapWave(int wid, int wave);
int IKUnmapWave(int wid, int wave);
int IKRtGrColors(int wid, int colors[]);
int IKAlGrColors(int wid, const int colors[]);
int IKSetOverlayColor(int wid, int icolor);
int IKAlSclBar(int wid, int flag);
int IKIsSclBarOn(int wid);
int IKRtSclBarAttr(
    int wid, int* color, int* direction, float* len, int* thick, int pos[]
);
int IKAlSclBarAttr(
    int wid, int color, int direction, float len, int thick, const int pos[]
);
int IKRtStepAttr(int wid, int attr[]);
int IKAlStepAttr(int wid, const int attr[]);
int IKRtDisOffset(int wid, int wave, int offset[]);
int IKAlDisOffset(int wid, int wave, const int offset[]);
int IKRtDisResLevel(int wid);
int IKAlDisResLevel(int wid, int istream);
int IKAlDisSec(int wid, int sec_pos, int check_mapped);
int IKRtDisSec(int wid);
int IKAlDisSec_ZWT(int wid, int sec, int wave, int time, int check_mapped);
int IKRtDisSec_ZWT(int wid, int* sec, int* wave, int* time);
char* IKLockSec(int wid, int sec);
int IKUnlockSec(int wid, int sec);
char* IKLockSec_ZWT(int wid, int sec, int wave, int time);
int IKUnlockSec_ZWT(int wid, int sec, int wave, int time);
int IKAlSecScl(int wid, int sec, const float scale[]);
int IKRtSecScl(int wid, int sec, float scale[]);
int IKRtSecMMM(int wid, int sec, float* min, float* max, float* mean);
int IKAlSecMMM(int wid, int sec, float min, float max, float mean);
int IKRtZWTMMM(
    int wid, int sec, int wave, int time, float* min, float* max, float* mean
);
int IKAlZWTMMM(
    int wid, int sec, int wave, int time, float min, float max, float mean
);
int IKAlScl(int wid, int wave, const float scale[]);
int IKRtScl(int wid, int wave, float scale[]);
int IKRtDisWave(int wid, int w_array[]);
int IKAlSWFile(int wid, const char* file_name);

int IKNextSec(int wid);
int IKPrevSec(int wid);
int IKLastSec(int wid);
int IKFirstSec(int wid);
int IKDisplay(int wid);
int IKRaiseWin(int wid);
int IKLowerWin(int wid);
int IKAllocTrueColor(int wid, int colors[]);
int IKClearWin(int wid);
int IKClearWave(int wid, int wave);
int IKRtFontAscent(int wid);
int IKWinToData(int wid, int wave, const int w_pos[], int d_pos[], int npts);
int IKDataToWin(int wid, int wave, int w_pos[], const int d_pos[], int npts);
int IKDataToDVPos(
    int wid, float ind_x, float ind_y, float ind_z, float rtn_pos[3]
);
int IKXEvtToData(
    int wid,int iwave, XEvent* event, float xyz[3], int* p_iw, int* p_it
);
int IKIncZoom(int wid);
int IKDecZoom(int wid);
int IKRtWaveGrColor(int wid, int wave);
int IKFirstImage(int wid, int wave, int check_mapped);
int IKLastImage(int wid, int wave, int check_mapped);
int IKNextImage(int wid, int wave, int inc, int check_mapped);
int IKRegMouseButton(int wid, Window send_event_win, int button_mask);
int IKUnRegMouseButton(int wid, Window send_event_win, int button_mask);
int IKRegEvt(int wid, Window send_event_win, long event_mask);
int IKUnregEvt(int wid, Window send_event_win, long event_mask);
int IKRegDisChg(int wid, Window send_event_win);
int IKUnregDisChg(int wid, Window send_event_win);
int IKInitEventTable(int* evtable);
int IKSignalMonitor(int wid, XEvent event);
int IKPopMemList(int wid);
int IKWinExistState(int wid);
int IKAttachWin(int istream, int wid, const char* status);
int IKAlFilename(int wid, int iwave, const char* filename);
int IKRtFilename(int wid, int iwave, char* filename);
int IKAlComplexDis(int wid, int flag);
int IKRtComplexDis(int wid);
int IKAlSclAlgorithm(int wid, int algorithm_code);
int IKRtSclAlgorithm(int wid);
int IKAlDisSyncMode(int wid, int mode);
int IKRtDisSyncMode(int wid);
int IKAlOffsetGroup(int wid, const int wave[]);
int IKRtOffsetGroup(int wid, int wave[]);
/* graphics routine */
int IKAlGrDis(int wid, int iflag);
int IKRtGrDis(int wid);
int IKGrAddPt(
    int wid,
    float x_coord,
    float y_coord,
    float z_coord,
    int wave,
    int time,
    int shape,
    int color,
    int type
);
int IKGrAddLn(
    int wid,
    float x1,
    float y1,
    float z1,
    float x2,
    float y2,
    float z2,
    int wave,
    int time,
    int thickness,
    int pattern,
    int color,
    int type
);
int IKGrAddBox(
    int wid,
    float x,
    float y,
    float z,
    int wave,
    int time, 
    float width,
    float height,
    int thickness,
    int pattern,
    int color,
    int type
);
int IKGrAddCir(
    int wid,
    float x_center,
    float y_center,
    float z_center,
    int wave,
    int time,
    float radius,
    int thickness,
    int color,
    int type
);
int IKGrAddStr(
    int wid,
    float x_start,
    float y_start,
    float z_start,
    int wave,
    int time,
    const char* in_string,
    int size,
    int color,
    int type
);
int IKGrAddMultStr(
    int wid,
    const IW_POINT_PTR plist,
    int wave,
    int time,
    int npts,
    int max_char,
    char* in_string[],
    int size,
    const int color[],
    int type
);
int IKGrAtt5DMod(int iwid);
int IKGrDet5DMod(int iwid);
int IKUpd5DMod(
    const float pt[],
    IW_GR_5D_MODEL_PTR mod_list,
    int* iwave,
    int itime,
    int flag
);
int IKGrGetPtr(int wid, int id, IW_GRL_PTR* ptr, int* sec);
int IKGrRmGrID(int wid,int id);
int IKGrRmGrList(int wid, const int* list, int nobj);
int IKGrRmProc(int wid);
int IKGrRmAllGrID(int wid);
int IKGrDisGrID(int wid, int igid);
int IKGrDisGrList(int wid, const int* ulist, int nobj);
int IKBgnFastGr(int wid);
int IKEndFastGr(int wid);
int IKDisFastGr(int wid);
int IKGrAddPts(
    int wid,
    const IW_POINT_PTR points_list,
    int num_points,
    int wave,
    int time,
    int shape, 
    int color,
    int type
);
int IKGrAddLns(
    int wid,
    const IW_POINT_PTR points_list,
    int num_points,
    int wave,
    int time,
    int thickness,
    int pattern,
    int color,
    int type
);
int IKGrAddMultLns(
    int wid,
    const IW_POINT_PTR points_list,
    int num_points,
    int wave,
    int time,
    int thickness,
    int pattern,
    int color,
    int type
);
int IKGrAddPoly(
    int wid,
    const IW_POINT_PTR points_list,
    int num_points,
    int wave,
    int time,
    int thickness,
    int color,
    int type
);
void IKDeleteGrNode(IVEPtrRep* p_head_off, IW_GRL_PTR node);
void IKInitImgMem(void);
int IKSaveImage(int wid, const char* f_name);
int IKCheckILImg(int wid, int sec, int res, int type);
void IKCompareRes(int wid, int sec, int res, int type);
int IKRmImgMem(int wid,int sec, int type);
int IKImgModified(int wid, int sec);
int IKBufModified(int wid, int sec, int res);
#ifdef _LARGEFILE64_SOURCE
off64_t IKGetResolutionOffset(int wid, int sec, int res, int use_pix);
#else
off_t IKGetResolutionOffset(int wid, int sec, int res, int use_pix);
#endif
int IKCreateWin(int wid);
int IKDeleteWin(int wid, int complete_delete, int new_is_scratch);
int IKClearOverlay(int wid);
int IKRtCurDataGraphics(int wid, int z_offset, IW_GRL_PTR** graphic_ptr);
int IKRtCurWinGraphics(int wid, IW_GRL_PTR** graphic_ptr);
int IKIsViewFile(int wid);
int IKRtStereoOffset(int wid);
int IKAlPointerPos(int wid, const int* xypos);
int IKAlCursor(int wid, int shape);
int IKInitFPSCounters(
    int wid, unsigned int* p_frame_counter, double* p_time_pt
);
int IKMeasureFPS(
    int wid,
    unsigned int frame_counter,
    double time_pt,
    double* p_fps,
    double* p_intvl_s
);
int IKParseWindowID(const char* name);
int IKRtThreshold(int wid, int wave, float* p_threshold);
int IKAlThreshold(int wid, int wave, float threshold);
int IKClrThreshold(int wid, int wave);
int IKRtSecThreshold(int wid, int section, float* p_threshold);
int IKAlSecThreshold(int wid, int section, float threshold);
int IKClrSecThreshold(int wid, int section);
int IKCaptureImage(
    int istream,
    int ix,
    int iy,
    int width,
    int height,
    int buffers,
    int format,
    int type,
    int scan_size,
    void* pixels
);
int IKAlSecNumColor(int wid, int color_index);
int IKRtSecNumColor(int wid);
int IKAlSecNumDisplay(int wid, int style);
int IKRtSecNumDisplay(int wid);

/*
 * To avoid extra copy operations, specially implement the calls to add
 * multiple graphics for Fortran.
 */
#ifdef IVE_HAS_FORTRAN
#include "ive_fortran.h"

int IKGrAddMultStrFtn(
    int wid,
    const ftn_real_t plist[],
    int wave,
    int time,
    int npts,
    int max_char,
    const ftn_char_t* labels,
    int size,
    const ftn_int_t colors[],
    int type,
    ftn_length_t labels_max_length
);
int IKGrAddPtsFtn(
    int wid,
    const ftn_real_t plist[],
    int num_points,
    int wave,
    int time,
    int shape,
    int color,
    int type
);
int IKGrAddLnsFtn(
    int wid,
    const ftn_real_t plist[],
    int num_points,
    int wave,
    int time,
    int thickness,
    int pattern,
    int color,
    int type
);
int IKGrAddMultLnsFtn(
    int wid,
    const ftn_real_t plist[],
    int num_points,
    int wave,
    int time,
    int thickness,
    int pattern,
    int color,
    int type
);
int IKGrAddPolyFtn(
    int wid,
    const ftn_real_t plist[],
    int num_points,
    int wave,
    int time,
    int thickness,
    int color,
    int type
);
#endif /* IVE_HAS_FORTRAN */

/* This is from IWNewMemMgr.c.  It is used in MonClient.c */
void FreeBufImg(int wid);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* include guard */
