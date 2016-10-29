#ifndef __CUDA_SIRECON_H
#define __CUDA_SIRECON_H
#ifdef _WIN32
#define _USE_MATH_DEFINES
#endif
#include <cuda.h>
#include <cufft.h>
#include <cuComplex.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <stdexcept>

#ifndef __clang__
#include <omp.h>
#endif

#include <memory>
#include <cfloat>
#include <cassert>

#include <complex>

#ifndef NDEBUG
#ifndef _WIN32
#include <nvml.h>  //for memory query
static nvmlDevice_t nvmldevice;
static nvmlMemory_t memoryStruct;
#endif
#endif

#include "Buffer.h"
#include "GPUBuffer.h"
#include "CPUBuffer.h"
#include "PinnedCPUBuffer.h"
#include "gpuFunctions.h"

#ifdef __SIRECON_USE_TIFF__
#include <tiffio.h>
#define cimg_use_tiff
#include <CImg.h>
using namespace cimg_library;
#else
#include <IMInclude.h>  // MRC file I/O routines
#endif

// Block sizes for reduction kernels
#define RED_BLOCK_SIZE_X 64
#define RED_BLOCK_SIZE_Y 4

// Some simulation parameters
/**ratio of beam size in pupil to pupil size */
#define SPOTRATIO 0.1
#define MAXPHASES 25
/** If k0 initial estimate is off from the guess by over this many pixels, a warning will be displayed */
#define K0_WARNING_THRESH 2


#define CHECKED_DELETE(PTR) \
  if (PTR) {\
    delete PTR;\
    PTR = 0;\
  }
#define CHECKED_DELETE_ARR(PTR) \
  if (PTR) {\
    delete [] PTR;\
    PTR = 0;\
  }                         

#ifdef __SIRECON_USE_TIFF__
static TIFF *otf_tiff;  // TODO: get rid of this; use "CImg<> m_otf_tiff" member in SIM_reconstructor instead.
#else
static const int istream_no = 1;
static const int ostream_no = 2;
static const int otfstream_no = 3;
static const int aligned_stream_no = 10;
static const int separated_stream_no = 11;
static const int overlaps_stream_no = 12;

// static IW_MRC_HEADER header;
static IW_MRC_HEADER aligned_header;
static IW_MRC_HEADER sep_header;
static IW_MRC_HEADER overlaps_header;

struct myExtHeader {
  float timestamp;
  float phaseAbsDeg;
  float expDose;
  float xdrift;
  float ydrift;
};
#endif

static float maxval = -FLT_MAX;
static float minval = FLT_MAX;


struct vector {
  float x;
  float y;
};
struct vector3d {
  float x;
  float y;
  float z;
};

/** Overall image reconstruction parameters. */
struct ReconParams {
  float k0startangle, linespacing;
  float na, nimm;
  int   ndirs, nphases, norders_output;
  int norders;
  float *phaseSteps; /** user-specified non-ideal phase steps, one for each orientation */
  int   bTwolens;    /** whether to process I{^5}S dataset */
  int   bFastSIM;   /** fast SIM data is organized differently */

  /* algorithm related parameters */
  float zoomfact;
  int   z_zoom;
  int   nzPadTo;  /** pad zero sections to this number of z sections */
  float explodefact;
  int   bFilteroverlaps;
  int   recalcarrays; /** whether to calculate the overlaping regions between bands just once or always; used in fitk0andmodamps() */
  int   napodize;
  std::vector<float> forceamp;
  // float *k0angles;
  std::vector<float> k0angles;
  int   bSearchforvector;
  int   bUseTime0k0;   /** whether to use time 0's k0 fit for the rest of a time series */
  int   apodizeoutput;  /** 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands() */
  float apoGamma;
  int   bSuppress_singularities;  /** whether to dampen the OTF values near band centers; used in filterbands() */
  int   suppression_radius;   /** if suppress_singularities is 1, the range within which suppresion is applied; used in filterbands() */
  int   bDampenOrder0;  /** whether to dampen the OTF values near band centers; used in filterbands() */
  int   bFitallphases;  /** In NL SIM, whether to use fitted phase for all orders or to infer higher order phase from order 1 phase fil */
  int   do_rescale; /** fading correction method: 0-no correction; 1-with correction */
  int   equalizez;
  int   equalizet;
  int   bNoKz0;   /** if true, no kz=0 plane is used in modamp fit and assemblerealspace() */
  float wiener, wienerInr;
  int   bUseEstimatedWiener;

  /* OTF specific parameters */
  int   nxotf, nyotf, nzotf;
  float dkzotf, dkrotf;  /** OTF's pixel size in inverse mirons */
  int   bRadAvgOTF;   /** is radially-averaged OTF used? */

  /* drift correction and phase step correction related flags */
  int   bFixdrift;   /** whether nor not to correct drift between pattern directions */
  float drift_filter_fact; /** fraction of the NA; used as the cutoff frequency in high-pass filtering in drift estimation */

  /* Camera related parameters */
  float constbkgd;
  int bBgInExtHdr; /** In Andor EMCCD, background varies with each exposure, esp. in EM mode. Hidden-behind-aluminum-foil pixels can be used to estimate background of each exposure and stored in the extended header. When this option is true, the 3rd float of the extended header stores such estimated background values. */
  int   bUsecorr;    /** whether to use a camera flat-fielding (or correction) file */
  char  corrfiles[400];  /** name of the camera correction file if bUsecorr is 1 */
  float readoutNoiseVar;
  float electrons_per_bit;

  /* Debugging flags */
  int   bMakemodel;  /** whether to fake an ideal point source and obtain a recontruction of it (i.e., an effective PSF in some sense) */
  int   bSaveSeparated; /** whether to save separated bands and then quit before filterbands */
  char  fileSeparated[400];
  int   bSaveAlignedRaw; /** whether to save dirft-corrected raw images (within each direction) and then quit before filterbands */
  char  fileRawAligned[400];
  int   bSaveOverlaps; /** whether to save makeoverlaps() output into a file */
  char  fileOverlaps[400];

  char ifiles[400];
  char ofiles[400];
  char otffiles[400];
  int ifilein;
  int ofilein;
  int otffilein;
};
struct ImageParams {
  int nx;
  int ny;
  int nz;
  int nz0;
  short nwaves;
  short wave[5];
  short ntimes;
  unsigned short curTimeIdx;
  float dx;
  float dy;
  float dz;
  float inscale;
};
struct DriftParams {
  vector3d* driftlist;
  float* phaseList;
  vector* driftAcq;
  vector3d* drifts;
  float* phaseAbs;
  float* timestamps;
  float* timestamps_for_fitting;
  float* expDose;
  vector3d* drift_bt_dirs;
  DriftParams() : driftlist(0), phaseList(0), driftAcq(0), drifts(0),
  phaseAbs(0), timestamps(0), timestamps_for_fitting(0), expDose(0) {
  };
  ~DriftParams() {
    CHECKED_DELETE_ARR(driftlist);
    CHECKED_DELETE_ARR(phaseList);
    CHECKED_DELETE_ARR(driftAcq);
    CHECKED_DELETE_ARR(drifts);
    CHECKED_DELETE_ARR(phaseAbs);
    CHECKED_DELETE_ARR(timestamps);
    CHECKED_DELETE_ARR(timestamps_for_fitting);
    CHECKED_DELETE_ARR(expDose);
  };
};
struct ReconData {
  int sizeOTF;
  std::vector<GPUBuffer> otf;
  CPUBuffer background;
  CPUBuffer slope;
  float backgroundExtra;
  std::vector<std::vector<GPUBuffer> > savedBands;
  std::vector<float> sepMatrix;
  std::vector<float> noiseVarFactors;
  GPUBuffer overlap0;
  GPUBuffer overlap1;
  std::vector<vector> k0;
  std::vector<vector> k0_time0;
  std::vector<vector> k0guess;
  std::vector<std::vector<cuFloatComplex> > amp;
  std::vector<double> sum_dir0_phase0;
  GPUBuffer bigbuffer;
  GPUBuffer outbuffer;
};


void SetDefaultParams(ReconParams *pParams);


// Functions for dealing with input images

/** General setup applicable to all times and waves
 * Sets up image parameters and sets up OTFs.
 * Allocates memory for images on GPU, and for sepmatrix and noisevars on CPU.
 * */
// void setup(ReconParams* params, ImageParams*
//     imgParams, DriftParams* driftParams, ReconData* data);
void setup_part2(ReconParams* params, ImageParams* imgParams, ReconData* reconData);
// #ifdef __SIRECON_USE_TIFF__
// void setup(CImg<> &inTIFF, ReconParams* params, ImageParams* imgParams, ReconData* reconData);
// #endif
void loadHeader(const ReconParams& params, ImageParams* imgParams, IW_MRC_HEADER &header);
// void readDriftData(const ReconParams& params, DriftParams* driftParams);
void getOTFs(ReconParams* params, const ImageParams& imgParams,
    ReconData* data);
void determine_otf_dimensions(int norders, int nz, ReconParams *params,
    int *sizeOTF);
void allocateOTFs(int norders, int sizeOTF,
    std::vector<GPUBuffer>* otfs);
int loadOTFs(int norders, const ReconParams& params,
    const ImageParams& imgParams, ReconData* data);
void allocateImageBuffers(const ReconParams& params,
    const ImageParams& imgParams, ReconData* reconData);

void setOutputHeader(const ReconParams& myParams, const ImageParams& imgParams,
                     IW_MRC_HEADER &header);

void bgAndSlope(const ReconParams& myParams,
    const ImageParams& imgParams, ReconData* reconData);

void getbg_and_slope(const char *corrfiles, float *background,
    float *slope, int nx, int ny);

void makematrix(int nphases, int norders, int dir, float *arrPhases,
    float *sepMatrix, float * noiseVarFactors);

void allocSepMatrixAndNoiseVarFactors(const ReconParams& params,
    ReconData* reconData);

#ifdef __SIRECON_USE_TIFF__
void load_and_flatfield(CImg<> &cimg, int section_no, float *bufDestiny, 
    float background, float inscale);
#endif

void load_and_flatfield(int section_no, int wave_no,
    int time_no, float *bufDestiny, float *buffer, int nx, int ny,
    float *background, float backgroundExtra, float *slope, float inscale,
    int bUsecorr);

void saveIntermediateDataForDebugging(const ReconParams& params);

void matrix_transpose(float* mat, int nRows, int nCols);

/*!
  
*/
void findModulationVectorsAndPhasesForAllDirections(
    int zoffset, ReconParams* params, const ImageParams& imgParams,
    DriftParams* driftParams, ReconData* data);

void apodizationDriver(int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data);
void rescaleDriver(int it, int iw, int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data);
void transformXYSlice(int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data);

int fitXYdrift(vector3d *drifts, float * timestamps, int nPoints,
    vector3d *fitted_drift, float *eval_timestamps, int nEvalPoints);
void calcPhaseList(float * phaseList, vector3d *driftlist,
    float *phaseAbs, float k0angle, float linespacing,
    float dr, int nphases, int nz, int direction, int z);

void writeResult(int it, int iw, const ReconParams& params,
    const ImageParams& imgParams, const ReconData& reconData);
#ifndef __SIRECON_USE_TIFF__
void saveCommandLineToHeader(int argc, char **argv, IW_MRC_HEADER &header);
#endif

void dumpBands(std::vector<GPUBuffer>* bands, int nx, int ny, int nz0);

void deviceMemoryUsage();

double meanAboveBackground_GPU(GPUBuffer &img, int nx, int ny, int nz);
void rescale_GPU(GPUBuffer &img, int nx, int ny, int nz, float scale);

// Compute rdistcutoff
int rdistcutoff(int iw, const ReconParams& params, const ImageParams& imgParams);
float get_phase(cuFloatComplex b);
float cmag(cuFloatComplex a);
cuFloatComplex cmul(cuFloatComplex a, cuFloatComplex b);

// Extern Blas and Lapack declarations
extern "C" void sgemm_(const char*, const char*, int*, int*, int*,
    float*, float*, int*, float*, int*, float*, float*, int*);
extern "C" void sgetrf_(int*, int*, float*, int*, int*, int*);
extern "C" void sgetri_(int*, float*, int*, int*, float*, int*, int*);
extern "C" void sgels_(const char*, int*, int*, int*, float*, int*, float*,
    int*, float*, int*, int*);

#ifdef __SIRECON_USE_TIFF__
extern "C" int load_tiff(TIFF *const tif, const unsigned int directory, const unsigned colind, float *const buffer);
extern "C" int save_tiff(TIFF *tif, const unsigned int directory, int colind, const int nwaves, int width, int height, float * buffer , int bIsComplex);

std::vector<std::string> gatherMatchingFiles(std::string target_path, std::string pattern);
std::string makeOutputFilePath(std::string inputFileName, std::string insert);
#endif

#endif
