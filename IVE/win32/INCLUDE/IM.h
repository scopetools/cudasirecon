#ifndef UCSF_MSG_IM_H
#define UCSF_MSG_IM_H

/* @(#) $Id: IM.h,v 1.2 2002/01/18 00:24:43 eric Exp $ */
/* $Name:  $ */
/*
 * Provides function prototypes for IM functions.
 */

#ifdef __cplusplus
extern "C" {
#endif

int IMOpen(int istream, const char* name, const char* attrib);
void IMClose(int istream);
int IMUnit(int istream);
void IMRdHdr(
    int istream,
    int ixyz[3],
    int mxyz[3],
    int* imode,
    float *min,
    float *max,
    float *mean
);
void IMWrHdr(
    int istream,
    const char title[80],
    int ntflag,
    float dmin,
    float dmax,
    float dmean
);
void IMGetExHdr(int istream, void* extended);
void IMGetHdr(int istream, void* header);
void IMPutHdr(int istream, const void* header);
void IMTrHdr(int istream, int jstream);
void IMTrCel(int istream, int jstream);
void IMTrExHdr(int jstream, int istream, int jsec, int isec);
int IMTrFmt(int istream, int jstream);
void IMTrLab(int istream, int jstream);
void IMCrHdr(
    int istream,
    const int ixyz[3],
    const int mxyz[3],
    int imode,
    const char* labels,
    int nl
);
void IMAlCel(int istream, const float cell[6]);
void IMAlCon(int istream, int flag);
void IMAlDat(
    int istream, int itype, int lensnum, int n1, int n2, float v1, float v2
);
void IMAlDel(int istream, const float delta[3]);
void IMAlDis(int flag);
void IMAlExHdr(int istream, int isec, const int ivals[], const float rvals[]);
void IMAlExHdrSize(int istream, int nints, int nreals, int nsecs);
void IMAlExHdrZWT(
    int istream, int iz, int iw, int it, const int ival[], const float rval[]
);
void IMAlExt(int istream, const void* extra, int istart, int nextra);
int IMAlFmt(int istream, int iformat);
void IMAlLab(int istream, const char* labels, int nl);
void IMAlMap(int istream, const int mcrs[3]);
void IMAlMode(int istream, int imode);
void IMAlOrig(int istream, float xorig, float yorig, float zorig);
void IMAlPrt(int flag);
void IMAlRes(int istream, int mres, int mzfact);
void IMAlSam(int istream, const int mxyz[3]);
void IMAlSiz(int istream, const int inxyz[3], const int nxyzst[3]);
void IMAlSpg(int istream, int nspg, int mbsym);
void IMAlSym(int istream, int mbsym, const void* jbsym);
void IMAlTlt(int istream, const float vals[3]);
void IMAlTltRot(int istream, const float vals[3]);
void IMAlTSt(int istream, int itst);
void IMAlWav(int istream, int nwave, const float wavelen[]);
void IMAlWavMM(int istream, int nwave, float dmin, float dmax);
void IMAlZWT(int istream,int nz, int nw, int nt, int wflag);
void IMRtCel(int istream, float cell[6]);
void IMRtDat(
    int istream,
    int* itype,
    int* lensnum,
    int* n1,
    int* n2,
    float* v1,
    float* v2
);
void IMRtDel(int istream, float delta[3]);
void IMRtExHdr(int istream, int jsec, int ivals[], float rvals[]);
void IMRtExHdrSize(int istream, int* nints, int* nreals);
void IMRtExHdrZWT(
    int istream, int iz, int iw, int it, int ival[], float rval[]
);
void IMRtExt(int istream, void* extra, int istart, int nextra);
void IMRtFmt(int istream, int* iformat);
void IMRtLab(int istream, char labels[20][80], int* nlabels);
void IMRtMap(int istream, int mcrs[3]);
void IMRtMode(int istream, int* mode);
void IMRtMst(int istream, int nxyzst[3]);
void IMRtOrig(int istream, float* xorig, float* yorig, float* zorig);
void IMRtRes(int istream, int* mres, int* mzfact);
void IMRtResInfo(
    int istream, int* nx, int* ny, int* nz, int* mxyfact, int* mzfact
);
void IMRtSam(int istream, int mxyz[3]);
void IMRtSiz(int istream, int inxyz[3], int mxyz[3], int nxyzst[3]);
void IMRtSpg(int istream, int* nspg, int* mbsym);
int IMRtStream(int wid);
void IMRtSym(int istream, int* mbsym, void* jbsym);
void IMRtTlt(int istream, float tilt[3]);
void IMRtTSt(int istream, int* itst);
void IMRtWav(int istream, int* nwave, float wavelen[]);
void IMRtWavMM(int istream, int iwave, float* dmin, float* dmax);
int IMRtWID(int istream);
void IMRtZWT(int istream, int* nz, int* nw, int* nt, int* wflag);
int IMPosn(int istream, int iz, int iy);
int IMPosnRes(int istream, int ires);
int IMPosnZWT(int istream, int iz, int iw, int it);
void IMRdLin(int istream, void* array);
void IMRdPal(int istream, void* array, const int* nx1, const int* nx2);
void IMRdPas(
    int istream,
    void* array,
    int mx,
    int my,
    int nx1,
    int nx2,
    int ny1,
    int ny2
);
void IMRdSec(int istream, void* array);
void IMRdSecl(int istream, void* array, const int* mlines);
void IMWrLin(int istream, const void* array);
void IMWrPal(int istream, const void* array, int nx1, int nx2);
void IMWrPas(
    int istream,
    const void* array,
    int mx,
    int my,
    int nx1,
    int nx2,
    int ny1,
    int ny2
);
void IMWrSec(int istream, const void* array);
void IMWrSecl(int istream, const void* array, int mlines);
void IMClLim(int istream, int xyzmin[3], int xyzmax[3], int mxyz[3]);
void IMFixExHdr(int istream, int izst, int izend, int izskip, int iflag);
void IMFixExHdrZWT(
    int istream,
    int izst,
    int izend,
    int nwt,
    const int iwtable[],
    int itst,
    int itend,
    int itinc,
    int iflag
);
void IMCalcDen(
    const float* array,
    int mx,
    int my,
    int nx1,
    int nx2,
    int ny1,
    int ny2,
    float* dmin,
    float* dmax,
    float* dmean
);
void IMInterp(
    const float* array_a,
    float* array_b,
    int nxa, 
    int nya,
    int nxb,
    int nyb,
    const float* amat,
    float xc,
    float yc,
    float xt,
    float yt,
    float scale
);

/***** new calls added by api *****/
int IMRtExHdrType(
    int iStreamNum, int* iExtHeaderType, int* iNumInts, int* iNumFloats
);
int IMRtExHdrValueZWT(
    int iStreamNum, int iZ, int iW, int iT, int iField, double* dValue
);
int IMAlExHdrValueZWT(
    int iStreamNum, int iZ, int iW, int iT, int iField, double dValue
);
int IMTrExHdr2(
    int iOutStream,
    int iInStream,
    int iz1,
    int iz2,
    int izinc,
    const int iProcWaveTable[IW_MAX_WAVE],
    int it1,
    int it2,
    int itinc
);
int IMRtZWTNum(int iStream, int* iz, int* iw, int* it, int iSecNum);
int IMRtSecNum(int iStream, int iz,  int iw, int it, int* iSecNum);


#ifdef __cplusplus
}
#endif

#endif /* include guard */
