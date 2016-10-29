#ifndef UCSF_MSG_IWL_H
#define UCSF_MSG_IWL_H

/*
 * The IW library provides all the IM function calls.
 */
#include "IM.h"
#include <X11/Xlib.h>


#ifdef __cplusplus
extern "C" {
#endif

int IWAttach(void);
int IWStart(int left, int bottom, int width, int height);
int IWStop(void);
const char* IWRetSHMDirectoryExt(void);
/* This is the deprecated version. */
int IWRetSHMDirectory(char* string_ptr);

int IWRtWID( int istream );
/* Display & Window Manipulation */
int IWRtDisArea(int istream,int *ix, int *iy, int *iwidth, int *iheight);
int IWRtWinPos(int istream, int *ixst, int *iyst, int *iwidth, int *iheight);
int IWAlWinPos(int istream, int ixst, int iyst, int iwidth, int iheight);
int IWRtMulDis(int istream, int *imdx, int *imdy);
int IWAlMulDis(int istream, int imdx, int imdy);
int IWRtDisFmt(int istream);
int IWAlDisFmt(int istream, int iorder);
int IWRtColorMode(int istream);
int IWAlColorMode(int istream, int mode);
int IWRtImgColor(int istream, int icolor[]);
int IWAlImgColor(int istream, const int icolor[]);
int IWRtGrDisRange(int istream, int range[]);
int IWAlGrDisRange(int istream, const int range[]);
int IWRtPsdGrColor(int istream);
int IWAlPsdGrColor(int istream, int ipsdgrcolor);
int IWAlDisWinBdr(int istream, int border_settings);
int IWRtDisWinBdr(int istream, int* p_border_settings);
int IWAlDisWinTool(int istream, int tool_settings);
int IWRtDisWinTool(int istream, int* p_tool_settings);
int IWRtTopImgLoc(int istream, int iwave);
int IWRtScratchSize(int istream);
int IWAlScratchSize(int istream, int isize);
int IWClrScr(int istream);
int IWAlClrBkg(int istream, int iflag);
int IWIsWinMap(int istream);
int IWMapWin(int istream);
int IWUnmapWin(int istream);
int IWRtDisDelay(int istream);
int IWAlDisDelay(int istream, int delay);
int IWRtInterp(int istream);
int IWAlInterp(int istream, int iflag);
int IWRtDisImg(int istream);
int IWAlDisImg(int istream, int iflag);
int IWRtBufImg(int istream);
int IWAlBufImg(int istream, int iflag);
int IWRtNoCopy(int istream);
int IWAlNoCopy(int istream, int iflag);
int IWAlFileFormat(int istream, int iflag);
int IWAttachFile(int istream, const char *file);
float IWRtZoom(int istream);
int IWAlZoom(int istream, float zoom);
int IWIsWaveMap(int istream, int iwave);
int IWAlWaveMap(int istream, int iwave, int iflag);
int IWRtGrColors(int istream, int icolor[]);
int IWAlGrColors(int istream, const int icolor[]);
int IWAlOverlay(int istream, int color);
int IWAlSclBar(int istream, int iflag);
int IWIsSclBarOn(int istream);
int IWRtSclBarAttr(int istream, int *icolor, int *idir, float *rlen,
                   int *ithick, int ipos[2]);
int IWAlSclBarAttr(int istream, int icolor, int idir,
                   float rlen, int ithick, const int ipos[2]);
int IWRtStepAttr(int istream, int iattr[2]);
int IWAlStepAttr(int istream, const int iattr[2]);
int IWRtDisOffset(int istream, int iwave, int ipixel[2]);
int IWAlDisOffset(int istream, int iwave, const int ipixel[2]);
int IWRtDisResLevel(int istream);
int IWAlDisResLevel(int istream, int ilevel);
int IWRtDisSec(int istream);
int IWAlDisSec(int istream, int isec);
int IWRtDisSec_ZWT(int istream, int *iz, int *iw, int *it);
int IWAlDisSec_ZWT(int istream, int iz, int iw, int it);
char *IWLockSec(int istream, int iz);
int IWUnlockSec(int istream, int iz);
char *IWLockSec_ZWT(int istream, int iz, int iw, int it);
int IWUnlockSec_ZWT(int istream, int iz, int iw, int it);
int IWRtSecScl(int istream, int isec, float scale[]);
int IWAlSecScl(int istream, int isec, const float scale[]);
int IWRtSecMMM(int istream, int isec, float *dmin, float *dmax, float *dmean);
int IWAlSecMMM(int istream, int isec, float dmin, float dmax, float dmean);
int IWRtZWTMMM(int istream, int iz, int iw, int it,
               float *dmin, float *dmax, float *dmean);
int IWAlZWTMMM(int istream, int iz, int iw, int it,
               float dmin, float dmax, float dmean);
int IWRtScl(int istream, int iwave, float scale[]);
int IWAlScl(int istream, int iwave, const float scale[]);
int IWRtDisWave(int istream, int *w_array);
int IWAlSWFile(int istream, const char *file);
int IWNextSec(int istream);
int IWPrevSec(int istream);
int IWLastSec(int istream);
int IWFirstSec(int istream);
int IWDisplay(int istream);
int IWRaiseWin(int istream);
int IWLowerWin(int istream);
int IWAllocTrueColor(int istream, int *icolor);
int IWClearWin(int istream);
int IWRtFontAscent(int istream);
int IWWinToData(int istream, int iwave, const int iwinpos[],
                int idatapos[], int npoint);
int IWDataToDVPos(int istream, float ind_x, float ind_y, float ind_z,
                                 float *rtn_pos);
int IWXEvtToDataExt(
    int istream, int iw, XEvent *event, float *xyz, int* p_iw, int *p_it
);
int IWXEvtToData(int istream, int iw, XEvent *event, float *xyz, int *p_it);
int IWIncZoom(int istream);
int IWDecZoom(int istream);
int IWRtWaveGrColor(int istream, int iwave);
int IWFirstImage(int istream, int iwave);
int IWLastImage(int istream, int iwave);
int IWNextImage(int istream, int iwave, int inc);
int IWPrevImage(int istream, int iwave, int inc);
int IWRegEvt(int istream, int ixwin, int mask);
int IWUnregEvt(int istream, int ixwin);
int IWRegDisChg(int istream, int ixwin);
int IWUnregDisChg(int istream, int ixwin);
int IWAlGrDis(int istream, int iflag);
int IWRtGrDis(int istream);
int IWGrAddPt(int istream, float x, float y, int shape, int icolor);
int IWGrAddPt2D(int istream, float x, float y, int shape, int icolor);
int IWGrAddPt3D(int istream, float x, float y, float z, int shape, int icolor);
int IWGrAddPt5D(int istream, float x, float y, float z,
                int wave, int time, int shape, int icolor);
int IWGrAddLn(int istream, float x1, float y1, float x2, float y2,
              int ith, int ipat, int icolor); 
int IWGrAddLn2D(int istream, float x1, float y1, float x2, float y2,
              int ith, int ipat, int icolor);
int IWGrAddLn3D(int istream, float x1, float y1, float z1, float x2, float y2,
              float z2, int ith, int ipat, int icolor);
int IWGrAddLn5D(int istream, float x1, float y1, float z1, float x2, float y2,
              float z2, int wave, int time, int ith, int ipat, int icolor);
int IWGrAddBox(int istream, float x, float y, float width, float height,
               int ith, int ipat, int icolor);
int IWGrAddBox2D(int istream, float x, float y, float width, float height,
               int ith, int ipat, int icolor);
int IWGrAddBox3D(int istream, float x, float y, float z,
                 float width, float height, int ith, int ipat, int icolor);
int IWGrAddBox5D(int istream, float x, float y, float z, int wave, int time,
                 float width, float height, int ith, int ipat, int icolor);
int IWGrAddCir(int istream, float cx, float cy, float rad, int ith, int icolor);
int IWGrAddCir2D(int istream, float cx, float cy, float rad, int ith, int icolor);
int IWGrAddCir3D(int istream, float cx, float cy, float cz, float rad, int ith, int icolor);
int IWGrAddCir5D(int istream, float cx, float cy, float cz,
                 int wave, int time, float rad, int ith, int icolor);
int IWGrAddStr(int istream, float xst, float yst, const char *str,
               int isiz, int icolor);
int IWGrAddStr2D(int istream, float xst, float yst, const char *str,
               int isiz, int icolor);
int IWGrAddStr3D(int istream, float xst, float yst, float zst, const char *str,
               int isiz, int icolor);
int IWGrAddStr5D(int istream, float xst, float yst, float zst,
                 int wave, int time, const char *str, int isiz, int icolor);
int IWGrAddMultStr(int istream, IW_POINT_2D_PTR plist, char *str[],
               int npts,int max_char, int isiz, const int icolor[]);
int IWGrAddMultStr2D(int istream, IW_POINT_2D_PTR plist, char *str[],
               int npts,int max_char, int isiz, const int icolor[]);
int IWGrAddMultStr3D(int istream, IW_POINT_PTR plist,
        char *str[], int npts,int max_char,int isiz, const int icolor[]);
int IWGrAddMultStr5D(int istream, IW_POINT_PTR plist, int wave,
                  int time, int npts, int max_char, char *str[],
                  int isiz, const int icolor[]);
int IWGrAtt5DMod(int istream);
int IWGrDet5DMod(int istream);
int IWUpd5DMod(float *pt, IW_GR_5D_MODEL_PTR mod_list, int *iwave,
              int itime, int flag);
int IWGrRmGrID(int istream, int igrid);
int IWGrRmGrList(int istream, const int list[], int nobj);
int IWGrRmProc(int istream);
int IWGrRmAllGrID(int istream);
int IWGrDisGrID(int istream, int igrid);
int IWGrDisGrList(int istream, const int list[], int nobj);
int IWBgnFastGr(int istream);
int IWEndFastGr(int istream);
int IWDisFastGr(int istream);
int IWGrAddPts(int istream, const IW_POINT_2D_PTR plist, int np,
                                              int shape, int icolor);
int IWGrAddPts2D(int istream, const IW_POINT_2D_PTR plist, int np,
                                              int shape, int icolor);
int IWGrAddPts3D(int istream, const IW_POINT_PTR plist, int np,
                                              int shape, int icolor);
int IWGrAddPts5D(int istream, const IW_POINT_PTR plist, int np,
                        int iwave, int itime, int shape, int icolor);
int IWGrAddLns(int istream, const IW_POINT_2D_PTR plist, int np, int ith,
                                              int ipat, int icolor);
int IWGrAddLns2D(int istream, const IW_POINT_2D_PTR plist, int np, int ith,
                                              int ipat, int icolor);
int IWGrAddLns3D(int istream, const IW_POINT_PTR plist, int np, int ith,
                                              int ipat, int icolor);
int IWGrAddLns5D(int istream, const IW_POINT_PTR plist, int np, int iwave,
                         int itime, int ith, int ipat, int icolor);
int IWGrAddMultLns(int istream, const IW_POINT_2D_PTR plist, int np,
                     int ith, int ipat, int icolor);
int IWGrAddMultLns2D(int istream, const IW_POINT_2D_PTR plist, int np,
                     int ith, int ipat, int icolor);
int IWGrAddMultLns3D(int istream, const IW_POINT_PTR plist, int np,
                     int ith, int ipat, int icolor);
int IWGrAddMultLns5D(int istream, const IW_POINT_PTR plist, int np,
                     int iwave, int itime, int ith, int ipat, int icolor);
int IWGrAddPoly(int istream, const IW_POINT_2D_PTR plist, int npt,
                  int ith, int icolor);
int IWGrAddPoly2D(int istream, const IW_POINT_2D_PTR plist, int npt,
                  int ith, int icolor);
int IWGrAddPoly3D(int istream, const IW_POINT_PTR plist, int npt,
                  int ith, int icolor);
int IWGrAddPoly5D(int istream, const IW_POINT_PTR plist, int npt,
                  int iwave, int itime, int ith, int icolor);
int IWWinExistState(int istream);
int IWAttachWin(int istream, int iwin, const char* status);
int IWDeleteWin(int istream, int mode);
int IWSaveImage(int istream, const char* f_name);
int IWDataToWin(
    int istream, int iwave, int iwinpos[], const int idatapos[], int npoint
);
int IWPopMemList(int istream);
int IWAlFilename(int istream, int iwave, const char* filename);
int IWRtFilename(int istream, int iwave, char* filename);
int IWRtFilenameLength(void);
int IWRtComplexDis(int istream);
int IWAlComplexDis(int istream, int iflag);
int IWAlPointerPos(int istream, const int winpos[2]);
int IWAlCursor(int istream,int shape);
int IWRtSclAlgorithm(int istream);
int IWAlSclAlgorithm(int istream, int algorithmCode);
int IWRtDisSyncMode(int istream);
int IWAlDisSyncMode(int istream, int mode);
int IWRtOffsetGroup(int istream, int wave[]);
int IWAlOffsetGroup(int istream, const int wave[]);
int IWRtMaxWin(void);
int IWRtTopDLWin(void);
int IWRtLastDLWin(void);
int IWRtNumWin(void);
int IWRtNumDisWin(void);
int IWRtScrnGeom(int* ixst, int* iyst, int* width, int* height);
int IWIsWinExist(int wid);
int IWRtNextDLWin(int wid);
int IWRtAllWin(int * nwin, int wlist[]);
int IWRtDisWinList(int* nwin, int wlist[]);
int IWRtColorNames(char* colors[IW_NUM_COLORS]);
int IWAlColorNames(char* colors[IW_NUM_COLORS]);
int IWIsIWLEvt( XEvent* event_ptr);
int IWGetButPressInfo(
    XEvent* event_ptr, int* wid_ptr, float* pxyz, int* pw, int* pt
);
int IWTransformPartialDataType(
    void* source,
    int width,
    int height,
    int mode,
    int x_st,
    int new_width,
    int y_st,
    int new_height,
    int new_mode,
    void* destination,
    int destination_stride
);
int IWCopyPartialData(
    void* source,
    int width,
    int height,
    int mode,
    int x_st,
    int new_width,
    int y_st,
    int new_height,
    void* destination,
    int destination_stride
);
int IWXEvtToWid(XEvent* event );
int IWXEvtToStream(XEvent* event );
int IWXEvtToKeySym(XEvent* event );
int IWRtLutSize(void);
void IWRtSclMinMax(int* imin, int* imax);
void IWRtLut(int ired[], int igreen[], int iblue[]);
void IWAlLut(const int ired[], const int igreen[], const int iblue[]);
int IWRtNumOVColors(void);
int IWAlOVColors(const int ired[], const int igreen[], const int iblue[]);
int IWRtOVColors(int ired[], int igreen[], int iblue[]);
int IWClearOverlay(int istream);
int IWStartOverlay(int istream, int icolor);
int IWEndOverlay(int istream);
int IWInvertMatrix(IW_MATRIX_PTR matrix,IW_MATRIX_PTR inv_matrix);
int IWVectorMatrixMult(
    IW_MATRIX matrix, IW_VECTOR vector, IW_VECTOR rtn_vector
);
int IWVectorVectorAdd(
    IW_VECTOR vector1, IW_VECTOR vector2, IW_VECTOR rtn_vector
);
int IWVectorVectorSub(
    IW_VECTOR vector1, IW_VECTOR vector2, IW_VECTOR rtn_vector
);
void IWMatrixFromAngles(float *angle,float matrix[3][3]);
void IWAnglesFromMatrix(float angle[3], float matrix[3][3]);
int IWScaleImage(
    const void* source_addr,
    const int xy[2],
    int mode,
    const float scl[4],
    IW_BYTE_MEM_PTR dest_addr
);
int IWScaleImageExt(
    const void* source_addr,
    const int xy[2],
    int mode,
    const float scl[4],
    int rmin,
    int rmax,
    IW_BYTE_MEM_PTR dest_addr
);
int IWRtPixSize(int mode);
void IWRtPixSizeExt(int mode, int* p_num, int* p_den, int* p_pad);
int IWIsStructSame(void* str1_ptr, void* str2_ptr, int size);
int IWIsValidMode(int mode);
float IWRtComplexNormalValue(IW_COMPLEX_PTR complex_ptr);
int IWSetComplexNormalValue(
    IW_COMPLEX_PTR complex_ptr, void* value_ptr, int value_mode
);
Display* IWGetAppXDisplay(void);
void IWSetAppXDisplay(Display *disp_id);
void* IWAlcSHM( int num, int size );
void* IWRealcSHM(void* ptr, int size);
void IWFreeSHM(void* ptr);
void* IWAlcNamedSHM(int num, int size, const char *name);
void* IWRtNamedSHMAdd(const char* name);
void IWFreeNamedSHM(void* temp_ptr);
IWEncodedSHMPtr IWEncodeSHMPtr(void* ptr);
void* IWDecodeSHMPtr(IWEncodedSHMPtr encoded_ptr);
/* IVEArena IWRetArenaPtr(void); */    /* need to include ive_shm.h */
void IWUSleep(unsigned int usecs);
const char* IWRetBaseDirectoryExt(void);
/* This is the deprecated version. */
int IWRetBaseDirectory(char *string_ptr);
const char* IWRetBinDirectory(void);
const char* IWRetVersion(void);
int IWDLFirstSec(void);
int IWDLLastSec(void);
int IWDLNextSec(void);
int IWDLPrevSec(void);
int IWDLIncZoom(void);
int IWDLDecZoom(void);
int IWFileExist(const char* filename);
int IWAlQuePri(int que0, int que1, int que2);
int IWRtCurDataGraphics(int istream, int z_offset, IW_GRL_PTR** graphic_ptr);
int IWRtCurWinGraphics(int istream, IW_GRL_PTR** graphic_ptr);
int IWIsViewFile(int istream);
int IWInitFPSCounters(
    int istream, unsigned int* p_frame_counter, double* p_time_pt
);
int IWMeasureFPS(
    int istream,
    unsigned int frame_counter,
    double time_pt,
    double* p_fps,
    double* p_intvl_s
);
int IWRtThreshold(int wid, int wave, float* p_threshold);
int IWAlThreshold(int wid, int wave, float threshold);
int IWClrThreshold(int wid, int wave);
int IWRtSecThreshold(int wid, int section, float* p_threshold);
int IWAlSecThreshold(int wid, int section, float threshold);
int IWClrSecThreshold(int wid, int section);
int IWCaptureImage(
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
int IWAlSecNumColor(int istream, int color_index);
int IWRtSecNumColor(int istream);
int IWAlSecNumDisplay(int istream, int style);
int IWRtSecNumDisplay(int istream);

#ifdef __cplusplus
}
#endif

#endif /* include guard */
