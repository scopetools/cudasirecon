
#include "cudaSirecon.h"
#include "cudaSireconImpl.h"
#include "SIM_reconstructor.hpp"

void SetDefaultParams(ReconParams *pParams)
{
  pParams->k0startangle =1.57193;
  pParams->linespacing = 0.177;  /* default to Nikon TIRF 100x */
  pParams->na=1.36;
  pParams->nimm=1.515;
  pParams->ndirs = 3;
  pParams->nphases = 3;
  pParams->phaseSteps = 0;
  pParams->norders_output = 0;
  pParams->bTwolens = 0;
  pParams->bFastSIM = 0;

  pParams->zoomfact = 2;
  pParams->z_zoom = 1;
  pParams->nzPadTo = 0;
  pParams->explodefact=1.0;
  pParams->bFilteroverlaps=1;
  pParams->recalcarrays = 1; /* whether to calculate the overlaping regions between bands just once or always; used in fitk0andmodamps() */
  pParams->napodize = 10;
  pParams->forceamp.assign(1, 0.f);
  // pParams->k0angles = NULL;
  pParams->bSearchforvector = 1;
  pParams->bUseTime0k0 = 1;  /* default to use time 0's k0 fit for the rest in a time series data */
  pParams->apodizeoutput = 0;  /* 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands() */
  pParams->apoGamma = 1;  /* 0-no apodize; 1-cosine apodize; 2-triangle apodize; used in filterbands() */
  pParams->bSuppress_singularities = 1;  /* whether to dampen the OTF values near band centers; used in filterbands() */
  pParams->suppression_radius = 10;   /* if suppress_singularities is 1, the range within which suppresion is applied; used in filterbands() */
  pParams->bDampenOrder0 = 0;  /* whether to dampen order 0 contribution; used in filterbands() */
  pParams->bFitallphases = 1;
  pParams->do_rescale=1; /* fading correction method: 0-no correction; 1-with correction */
  pParams->equalizez = 0;
  pParams->equalizet = 0;
  pParams->bNoKz0 = 0;
  pParams->bUseEstimatedWiener = 1;

  pParams->bRadAvgOTF = 0;  /* default to use non-radially averaged OTFs */

  pParams->bFixdrift = 0;
  pParams->drift_filter_fact = 0.0;

  pParams->constbkgd = 0.0;
  pParams->bBgInExtHdr = 0;
  pParams->bUsecorr = 0;
  pParams->corrfiles[0] = '\0';  /* name of CCD correction file if usecorr is 1 */
  pParams->electrons_per_bit = 0.6528;
  pParams->readoutNoiseVar = 32.42;  // electron^2

  pParams->bMakemodel = 0;
  pParams->bSaveAlignedRaw = 0;
  pParams->fileSeparated[0] = '\0';
  pParams->bSaveSeparated = 0;
  pParams->fileRawAligned[0] = '\0';
  pParams->bSaveOverlaps = 0;
  pParams->fileOverlaps[0] = '\0';

  pParams->ifilein = 0;
  pParams->ofilein = 0;
  pParams->otffilein = 0;
}


// Functions for dealing with images

bool notGoodDimension(unsigned num)
// Good dimension is defined as one that can be fatorized into 2s, 3s, 5s, and 7s
// According to CUFFT manual, such dimension would warranty fast FFT
{
  if (num==2 || num==3 || num==5 || num==7)
    return false;
  else if (num % 2 == 0) return notGoodDimension(num / 2);
  else if (num % 3 == 0) return notGoodDimension(num / 3);
  else if (num % 5 == 0) return notGoodDimension(num / 5);
  else if (num % 7 == 0) return notGoodDimension(num / 7);
  else
    return true;
}

unsigned findOptimalDimension(unsigned inSize, int step=-1)
{
  unsigned outSize = inSize;
  while (notGoodDimension(outSize))
    outSize += step;

  return outSize;
}

void setup_part2(ReconParams* params, ImageParams* imgParams, ReconData* reconData)
{
#ifdef __SIRECON_USE_TIFF__
  imgParams->nz0 = findOptimalDimension(imgParams->nz);
#endif

  getOTFs(params, *imgParams, reconData);
  allocSepMatrixAndNoiseVarFactors(*params, reconData);
  makematrix(params->nphases, params->norders, 0, 0,
      &(reconData->sepMatrix[0]), &(reconData->noiseVarFactors[0]));

  allocateImageBuffers(*params, *imgParams, reconData);

  imgParams->inscale = 1.0 / (imgParams->nx * imgParams->ny * imgParams->nz0 *
      params->zoomfact * params->zoomfact * params->z_zoom * params->ndirs);
  reconData->k0 = std::vector<vector>(params->ndirs);
  reconData->k0_time0 = std::vector<vector>(params->ndirs);
  reconData->k0guess = std::vector<vector>(params->ndirs);
  float delta_angle = M_PI / params->ndirs;
  float dkr = 1 / (imgParams->ny * imgParams->dy);  // assuming square images sizes
  float k0magguess = (1.0 / params->linespacing) / dkr;
  if (imgParams->nz > 1) {
    int nordersIn = params->nphases / 2 + 1;
    k0magguess /= nordersIn - 1; 
  }
  for (int i = 0; i < params->ndirs; ++i) {
    float k0angleguess;
    if (params->k0angles.size() < params->ndirs) {
      k0angleguess = params->k0startangle + i * delta_angle;
    } else {
      k0angleguess = params->k0angles[i];
    }
    reconData->k0guess[i].x = k0magguess * cos(k0angleguess);
    reconData->k0guess[i].y = k0magguess * sin(k0angleguess);
  }

  reconData->sum_dir0_phase0 = std::vector<double>(imgParams->nz *
      imgParams->nwaves);
  reconData->amp = std::vector<std::vector<cuFloatComplex> >(
      params->ndirs, std::vector<cuFloatComplex>(params->norders));
  for (int i = 0; i < params->ndirs; ++i) {
    reconData->amp[i][0].x = 1.0f;
    reconData->amp[i][0].y = 0.0f;
  }
}

#ifndef __SIRECON_USE_TIFF__
void loadHeader(const ReconParams& params, ImageParams* imgParams, IW_MRC_HEADER &header)
{
  int ixyz[3];
  int mxyz[3];
  int pixeltype;
  float min;
  float max;
  float mean;
  IMRdHdr(istream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(istream_no, &header);
  imgParams->nx = header.nx;
  imgParams->ny = header.ny;
  imgParams->nz = header.nz / (header.num_waves * header.num_times);
  /* header.nz is defined as the total number of sections =
   * nz*nwaves*ntimes (* ndirs*nphases in our case) */
  imgParams->nz /= params.nphases  * params.ndirs;
  imgParams->nwaves = header.num_waves;
  imgParams->ntimes = header.num_times;
  imgParams->wave[0]=header.iwav1;
  imgParams->wave[1]=header.iwav2;
  imgParams->wave[2]=header.iwav3;
  imgParams->wave[3]=header.iwav4;
  imgParams->wave[4]=header.iwav5;
  /* dy: lateral pixel size; dz: axial pixel size; both in microns */
  imgParams->dy = header.ylen;
  imgParams->dz = header.zlen;

  /* Initialize headers for intermediate output files if requested */
  if (params.bSaveAlignedRaw) {
    memcpy(&aligned_header, &header, sizeof(header));
    IMOpen(aligned_stream_no, params.fileRawAligned, "new");
    aligned_header.mode = IW_FLOAT;
    aligned_header.inbsym = 0;
    IMPutHdr(aligned_stream_no, &aligned_header);
  }
  if (params.bSaveSeparated) {
    memcpy(&sep_header, &header, sizeof(header));
    IMOpen(separated_stream_no, params.fileSeparated, "new");
    sep_header.nx = (imgParams->nx+2)/2;    // saved will be separated FFTs
    sep_header.mode = IW_COMPLEX;
    sep_header.inbsym = 0;
    IMPutHdr(separated_stream_no, &sep_header);
  }
  if (params.bSaveOverlaps) {
    memcpy(&overlaps_header, &header, sizeof(header));
    IMOpen(overlaps_stream_no, params.fileOverlaps, "new");
    overlaps_header.nz = imgParams->nz*2*params.ndirs*imgParams->ntimes*imgParams->nwaves;
    overlaps_header.num_waves = 2;  // save overlap 0 and 1 as wave 0 and 1 respectively
    overlaps_header.interleaved = WZT_SEQUENCE;
    overlaps_header.mode = IW_COMPLEX;   // saved will be full-complex overlaps in real space
    overlaps_header.inbsym = 0;
    IMPutHdr(overlaps_stream_no, &overlaps_header);
  }
  if (params.nzPadTo) {
    imgParams->nz0 = params.nzPadTo;
  } else {
    imgParams->nz0 = imgParams->nz;
  }
  printf("nx=%d, ny=%d, nz=%d, nz0 = %d, nwaves=%d, ntimes=%d\n",
      imgParams->nx, imgParams->ny, imgParams->nz, imgParams->nz0,
      imgParams->nwaves, imgParams->ntimes);
}
#endif

void getOTFs(ReconParams* params, const ImageParams& imgParams,
    ReconData* data)
{
  params->norders = 0;
  if (params->norders_output != 0) {
    params->norders = params->norders_output;
  } else {
    params->norders = params->nphases / 2 + 1;
  }
  determine_otf_dimensions(params->norders, imgParams.nz, params, &(data->sizeOTF));
  allocateOTFs(params->norders, data->sizeOTF, &(data->otf));
 loadOTFs(params->norders, *params, imgParams, data);
}

void determine_otf_dimensions(int norders, int nz, ReconParams *pParams, int *sizeOTF)
{
#ifdef __SIRECON_USE_TIFF__
  uint32 nxotf, nyotf, nzotf;
  float xres, yres;
  TIFFGetField(otf_tiff, TIFFTAG_IMAGEWIDTH, &nxotf);
  TIFFGetField(otf_tiff, TIFFTAG_IMAGELENGTH, &nyotf);
  TIFFGetField(otf_tiff, TIFFTAG_XRESOLUTION, &xres);
  TIFFGetField(otf_tiff, TIFFTAG_YRESOLUTION, &yres);
  nzotf = 0;
  do ++nzotf; while (TIFFReadDirectory(otf_tiff));

  /* determine nzotf, nxotf, nyotf, dkrotf, dkzotf based on dataset being 2D/3D and 
     flag bRadAvgOTF */

  if (nz == 1) {  /* 2D */
    pParams->nxotf = nxotf;
    if (pParams->bRadAvgOTF)
      pParams->nyotf = 1;
    else
      pParams->nyotf = nyotf;
    pParams->nzotf = 1;
    pParams->dkrotf = xres;  // dkrotf's unit is 1/micron
    pParams->dkzotf = 1;
  }
  else {   /* 3D */
    if (pParams->bRadAvgOTF) {
      pParams->nzotf = nxotf;
      pParams->nxotf = nyotf;
      pParams->nyotf = 1;
      pParams->dkzotf = xres;
      pParams->dkrotf = yres;
    }
    else {
      pParams->nzotf = nzotf / norders; // each order has a 3D OTF stack (non-negative kx half of Fourier space)
      pParams->nxotf = nxotf;
      pParams->nyotf = nyotf;
      pParams->dkzotf = xres;
      pParams->dkrotf = yres;
    }
  }
#else
  int ixyz[3], mxyz[3], pixeltype;
  float min, max, mean;
  IW_MRC_HEADER otfheader;
  /* Retrieve OTF file header info */
  IMRdHdr(otfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(otfstream_no, &otfheader);
  IMAlCon(otfstream_no, 0);
  /* determine nzotf, nxotf, nyotf, dkrotf, dkzotf based on dataset
   * being 2D/3D and flag bRadAvgOTF */
  if (nz == 1) {  /* 2D */
    pParams->nxotf = otfheader.nx;
    if (pParams->bRadAvgOTF)
      pParams->nyotf = 1;
    else
      pParams->nyotf = otfheader.ny;
    pParams->nzotf = 1;
    pParams->dkrotf = otfheader.xlen;  // dkrotf's unit is 1/micron
    pParams->dkzotf = 1;
  } else {   /* 3D */
    if (pParams->bRadAvgOTF) {
      pParams->nzotf = otfheader.nx;
      pParams->nxotf = otfheader.ny;
      pParams->nyotf = 1;
      pParams->dkzotf = otfheader.xlen;
      pParams->dkrotf = otfheader.ylen;
    } else {
      // each order has a 3D OTF stack (non-negative kx half of Fourier
      // space)
      pParams->nzotf = otfheader.nz / norders;
      pParams->nxotf = otfheader.nx;
      pParams->nyotf = otfheader.ny;
      pParams->dkzotf = otfheader.zlen;
      pParams->dkrotf = otfheader.xlen;
    }
  }
#endif
  printf("nzotf=%d, dkzotf=%f, nxotf=%d, nyotf=%d, dkrotf=%f\n",
      pParams->nzotf, pParams->dkzotf, pParams->nxotf, pParams->nyotf,
      pParams->dkrotf);

  /* sizeOTF and norders are determined so that correct memory can be
   * allocated for otf */
  *sizeOTF = pParams->nzotf*pParams->nxotf*pParams->nyotf;
}

void allocateOTFs(int norders, int sizeOTF,
    std::vector<GPUBuffer>* otfs)
{
  otfs->clear();
  for (int i = 0; i < norders; ++i) {
    GPUBuffer buff;
    buff.resize(sizeOTF * sizeof(cuFloatComplex));
    otfs->push_back(buff);
  }
}

int loadOTFs(int norders, const ReconParams& params,
    const ImageParams& imgParams, ReconData* data)
{
#ifdef __SIRECON_USE_TIFF__
  int i, nzotf;
  float *realpart, *imagpart;
  nzotf = 0;
  TIFFSetDirectory(otf_tiff, 0);
  do ++nzotf; while (TIFFReadDirectory(otf_tiff));

  realpart = (float *) malloc(data->sizeOTF * sizeof(float));
  imagpart = (float *) malloc(data->sizeOTF * sizeof(float));

  cuFloatComplex * otfTmp= (cuFloatComplex *) malloc(data->sizeOTF * sizeof(cuFloatComplex));
  CPUBuffer otfTmpBuffer(data->sizeOTF * sizeof(cuFloatComplex));

  /* Load OTF data, no matter 2D, 3D, radially averaged or not. */
  for (i=0; i<norders; i++) {
    int j;
    /* If OTF file has multiple sections, then read them into otf[i]; */
    if (imgParams.nz == 1 || params.bRadAvgOTF) {
      if (nzotf > i) { /* each section in OTF file is OTF of one order; so load that section into otf[i]  */
        load_tiff(otf_tiff, i, 0, realpart);
        load_tiff(otf_tiff, i, 1, imagpart);
        for (j=0; j<data->sizeOTF; j++) {
          otfTmp[j].x = realpart[j];
          otfTmp[j].y = imagpart[j];
        }
        otfTmpBuffer.setFrom((void*) otfTmp, 0, data->sizeOTF * sizeof(cuFloatComplex), 0);
        data->otf[i].setFrom(otfTmpBuffer, 0, data->sizeOTF * sizeof(cuFloatComplex), 0);
      }
      else   /* If there's just 1 OTF image, do not read any more and just duplicate otf[0] into otf[i] */
        data->otf[0].set(&(data->otf[i]), 0,
            data->otf[0].getSize(), 0);
    }
    else {  // non-radially averaged 3D OTF
      int nxyotf = params.nxotf * params.nyotf, z;
      for (z=0; z < params.nzotf; z++) {
        load_tiff(otf_tiff, z+i*nxyotf, 0, realpart);
        load_tiff(otf_tiff, z+i*nxyotf, 1, imagpart);
        for (j=0; j<nxyotf; j++) {
          otfTmp[j].x = realpart[j];
          otfTmp[j].y = imagpart[j];
        }
        otfTmpBuffer.setFrom((void*) otfTmp, 0, nxyotf * sizeof(cuFloatComplex), 0);
        data->otf[i].setFrom(otfTmpBuffer, 0, nxyotf * sizeof(cuFloatComplex),
                             z * nxyotf * sizeof(cuFloatComplex));
      }
    }
  }
  TIFFClose(otf_tiff);
  free(realpart);
  free(imagpart);
#else
  int ixyz[3], mxyz[3], pixeltype;
  float min, max, mean;
  IW_MRC_HEADER otfheader;
  /* Retrieve OTF file header info */
  IMRdHdr(otfstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(otfstream_no, &otfheader);

  CPUBuffer otfTmp(data->sizeOTF * sizeof(cuFloatComplex));

  /* Load OTF data, no matter 2D, 3D, radially averaged or not. */
  for (int i = 0; i < norders; i++) {
    /* If OTF file has multiple sections, then read them into otf[i]; */
    if (imgParams.nz == 1 || params.bRadAvgOTF) {
      if (otfheader.nz > i) {
        /* each section in OTF file is OTF of one order; so load that
         * section into otf[i]  */
        IMRdSec(otfstream_no, otfTmp.getPtr());
        otfTmp.set(&(data->otf[i]), 0, otfTmp.getSize(), 0);
      } else {
        /* If there's just 1 OTF image, do not read any more and just
         * duplicate otf[0] into otf[i] */
        data->otf[0].set(&(data->otf[i]), 0,
            data->otf[0].getSize(), 0);
      }
    } else {  // non-radially averaged 3D OTF
      for (int z = 0; z < params.nzotf; ++z) {
        IMRdSec(otfstream_no, otfTmp.getPtr());
        otfTmp.set(&(data->otf[i]),
            0, params.nxotf * params.nyotf * sizeof(cuFloatComplex),
            z * params.nxotf * params.nyotf * sizeof(cuFloatComplex));
      }
    }
  }
  IMClose(otfstream_no);
#endif

#ifndef NDEBUG
  for (std::vector<GPUBuffer>::iterator i = data->otf.begin();
      i != data->otf.end(); ++i) {
    assert(i->hasNaNs() == false);
  }
#endif
  return 1;
}

void allocateImageBuffers(const ReconParams& params,
    const ImageParams& imgParams, ReconData* data)
{
  data->savedBands.clear();
  data->savedBands.reserve(params.ndirs);
  for (int i = 0; i < params.ndirs; ++i) {
    data->savedBands.push_back(std::vector<GPUBuffer>());
    data->savedBands[i].reserve(params.nphases);
    for (int j = 0; j < params.nphases; ++j) {
      data->savedBands[i].push_back(GPUBuffer(
            (imgParams.nx / 2 + 1) * imgParams.ny * imgParams.nz0 *
            sizeof(cuFloatComplex), 0));
    }
  }
  for (int i = 0; i < params.ndirs; ++i) {
    for (int j = 0; j < params.nphases; ++j) {
      data->savedBands[i][j].setToZero();
    }
  }
}

#ifndef __SIRECON_USE_TIFF__
void setOutputHeader(const ReconParams& myParams, const ImageParams& imgParams,
                     IW_MRC_HEADER &header)
{
  header.mode = IW_FLOAT;
  header.nz = imgParams.nz * imgParams.nwaves * imgParams.ntimes *
    myParams.z_zoom;
  header.nx *= myParams.zoomfact;
  header.ny *= myParams.zoomfact;
  header.xlen /= myParams.zoomfact;
  header.ylen /= myParams.zoomfact;
  header.zlen /= myParams.z_zoom;
  header.inbsym = 0;
  IMPutHdr(ostream_no, &header);
  IMAlCon(ostream_no, 0);
}
#endif

void bgAndSlope(const ReconParams& myParams,
    const ImageParams& imgParams, ReconData* reconData)
{
  reconData->background.resize(sizeof(float) * imgParams.nx *
      imgParams.ny);
  reconData->slope.resize(sizeof(float) * imgParams.nx *
      imgParams.ny);
#ifndef __SIRECON_USE_TIFF__
  if (myParams.bUsecorr) {
    // flatfield correction of measured data using calibration data
    printf("loading CCD calibration file\n");
    getbg_and_slope(myParams.corrfiles,
        (float*)reconData->background.getPtr(),
        (float*)reconData->slope.getPtr(), imgParams.nx, imgParams.ny);
  } else {
#endif
    for (int i = 0; i < imgParams.nx * imgParams.ny; i++) {
      /* use the constant background value given by user */
      ((float*)(reconData->background.getPtr()))[i] = myParams.constbkgd;
      ((float*)(reconData->slope.getPtr()))[i] = 1.0;
    }
#ifndef __SIRECON_USE_TIFF__
  }
#endif
  reconData->backgroundExtra = 0;
}

#ifndef __SIRECON_USE_TIFF__
void getbg_and_slope(const char *corrfiles, float *background,
    float *slope, int nx, int ny)
{
  int cstream_no=10;
  int ixyz[3], mxyz[3], pixeltype;      /* variables for IMRdHdr call */
  float min, max, mean;      /* variables for IMRdHdr call */
  IW_MRC_HEADER header;

  if (IMOpen(cstream_no, corrfiles, "ro")) {
    fprintf(stderr, "File %s does not exist\n", corrfiles);
    exit(-1);
  }

  IMRdHdr(cstream_no, ixyz, mxyz, &pixeltype, &min, &max, &mean);
  IMGetHdr(cstream_no, &header);
  if (header.nx != nx || header.ny != ny) {
    fprintf(stderr,
        "calibration file %s has different dimension than data file",
        corrfiles);
    exit(-1);
  }
  IMAlCon(cstream_no, 0);

  IMRdSec(cstream_no, background);
  IMRdSec(cstream_no, slope);

  IMClose(cstream_no);
}
#endif

void findModulationVectorsAndPhasesForAllDirections(
    int zoffset, ReconParams* params, const ImageParams& imgParams,
    DriftParams* driftParams, ReconData* data)
{
  // Apodize (or edge softening) every 2D slice:
  apodizationDriver(zoffset, params, imgParams, driftParams, data);

  // 2D FFT every 2D slice:
  transformXYSlice(zoffset, params, imgParams, driftParams, data);

  for (int direction = 0; direction < params->ndirs; ++direction) {

    std::vector<GPUBuffer>* rawImages = &(data->savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(data->savedBands[direction]);

    std::vector<float> phaseList(params->nphases);
    if (params->phaseSteps != 0) {
      // User specified the non-ideal (i.e., not 2pi/5) phase steps for
      // each orientation. Under this circumstance, calculate the
      // non-ideal sepMatrix:
      for (int i = 0; i < params->nphases; ++i) {
        phaseList[i] = i * params->phaseSteps[direction];
      }
      makematrix(params->nphases, params->norders, direction,
          &(phaseList[0]), (&data->sepMatrix[0]),
          &(data->noiseVarFactors[0]));
    }

    // Unmixing info components in real or reciprocal space:
    separate(imgParams.nx, imgParams.ny, imgParams.nz,
        direction, params->nphases, params->norders,
        rawImages, &data->sepMatrix[0]);

#ifndef NDEBUG
    for (int phase = 0; phase < params->nphases; ++phase) {
      assert(rawImages->at(phase).hasNaNs(true) == false);
      std::cout << "Phase " << phase << "ok." << std::endl;;
    }
#endif

    if (imgParams.nz > 1) {
      // 1D FFT of a stack of 2D FFTs to obtain equivalent of 3D FFT
      cufftHandle fftplan1D;
      int fftN[1] = {imgParams.nz0};
      int stride = (imgParams.nx / 2 + 1) * imgParams.ny;
      int dist = 1;
      cufftResult cuFFTErr = cufftPlanMany(&fftplan1D, 1, fftN,
                                           fftN, stride, dist,
                                           fftN, stride, dist,
                                           CUFFT_C2C,
                                           (imgParams.nx / 2 + 1) * imgParams.ny);
      if (cuFFTErr != CUFFT_SUCCESS) {
        throw std::runtime_error("CUFFT plan creation failed");
      }
      for (int i = 0; i < params->nphases; ++i) {
        cufftExecC2C(fftplan1D, (cuFloatComplex*)((*bands)[i]).getPtr(),
                     (cuFloatComplex*)((*bands)[i]).getPtr(), CUFFT_FORWARD);
      }
      cufftDestroy(fftplan1D);
    }

    if (params->bMakemodel) {
      /* Use the OTF to simulate an ideal point source; replace bands
       * with simulated data
       * DM: k0 is not initialized but has memory allocate at this point.
       * k0 initialization code is near lines 430ff in sirecon.c */
      makemodeldata(imgParams.nx, imgParams.ny, imgParams.nz0, bands,
          params->norders, data->k0[direction], imgParams.dy, imgParams.dz,
          &data->otf, imgParams.wave[0], params);
    }

#ifndef __SIRECON_USE_TIFF__
    /* save the separated raw if requested */
    if (params->bSaveSeparated) {
      CPUBuffer tmp((*rawImages)[0].getSize());
      for (int phase = 0; phase < params->nphases; ++ phase) {
        (*rawImages)[phase].set(&tmp, 0, tmp.getSize(), 0);
        for (int z = 0; z < imgParams.nz; ++z) {
          float* imgPtr = (float*)tmp.getPtr();
          IMWrSec(separated_stream_no,
              imgPtr + (z + zoffset) * (imgParams.nx + 2) * imgParams.ny);
        }
      }
      continue; // skip the k0 search and modamp fitting
    }
#endif

    /* After separation and FFT, the std::vector rawImages, now referred to as
     * the std::vector bands, contains the center band (bands[0])
     * , the real (bands[1], bands[3], etc.) and
     * imaginary part (bands[2], bands[4], etc.) bands
     * of the different orders */

    if (! (imgParams.ntimes > 1 && imgParams.curTimeIdx > 0
           && params->bUseTime0k0) ) {
      data->k0[direction] = data->k0guess[direction];
      printf("k0guess[direction %d] = (%f, %f)\n", direction,
          data->k0guess[direction].x, data->k0guess[direction].y);
    }

    // Now to fix 3D drift between dirs estimated by determinedrift_3D()
    if (direction != 0 && params->bFixdrift) {
      fixdrift_bt_dirs(bands, params->norders, driftParams->drift_bt_dirs[direction],
          imgParams.nx, imgParams.ny, imgParams.nz0);
    }

    /* assume k0 vector not well known, so fit for it */
    cuFloatComplex amp_inv;
    cuFloatComplex amp_combo;
    if (params->bSearchforvector &&
        !(imgParams.ntimes > 1 && imgParams.curTimeIdx > 0 && params->bUseTime0k0)) {
      /* In time series, can choose to use the time 0 k0 fit for the
       * rest of the series.
       * Find initial estimate of modulation wave vector k0 by
       * cross-correlation. */
      findk0(bands, &data->overlap0, &data->overlap1, imgParams.nx,
          imgParams.ny, imgParams.nz0, params->norders,
          &(data->k0[direction]), imgParams.dy, imgParams.dz, &data->otf,
          imgParams.wave[0], params);

      if (params->bSaveOverlaps) {
        // output the overlaps
#ifdef __SIRECON_USE_TIFF__
        CPUBuffer tmp0(data->overlap0.getSize()*2);
        data->overlap0.set(&tmp0, 0, data->overlap0.getSize(), 0);
        data->overlap1.set(&tmp0, 0, data->overlap1.getSize(), data->overlap0.getSize());
        CImg<> ovlp0((float* )tmp0.getPtr(), imgParams.nx*2, imgParams.ny, imgParams.nz*2, 1,
          true);  // ovlp0 shares buffer with tmp0 (hence "true" in the final parameter)
        ovlp0.save_tiff(params->fileOverlaps);
#else
        CPUBuffer tmp0(data->overlap0.getSize());
        data->overlap0.set(&tmp0, 0, tmp0.getSize(), 0);
        CPUBuffer tmp1(data->overlap1.getSize());
        data->overlap1.set(&tmp1, 0, tmp1.getSize(), 0);
        cuFloatComplex* ol0Ptr = (cuFloatComplex*)tmp0.getPtr();
        cuFloatComplex* ol1Ptr = (cuFloatComplex*)tmp1.getPtr();
        for (int z = 0; z < imgParams.nz0; ++z) {
          IMWrSec(overlaps_stream_no, ol0Ptr + z * imgParams.nx * imgParams.ny);
          IMWrSec(overlaps_stream_no, ol1Ptr + z * imgParams.nx * imgParams.ny);
        }
#endif
      }
      printf("Initial guess by findk0() of k0[direction %d] = (%f,%f)\n", 
          direction, data->k0[direction].x, data->k0[direction].y);

      /* refine the cross-corr estimate of k0 by a search for best k0
       * vector direction and magnitude using real space waves*/

      printf("before fitk0andmodamp\n");
      /* if (myParams.recalcarrays==0 && dist>1.0) */
      params->recalcarrays = 1;
      /* if k0 is very close to the guess, we can save time by not
       * recalculating the overlap arrays */

      fitk0andmodamps(bands, &data->overlap0, &data->overlap1, imgParams.nx,
          imgParams.ny, imgParams.nz0, params->norders, &(data->k0[direction]),
          imgParams.dy, imgParams.dz, &data->otf, imgParams.wave[0],
          &data->amp[direction][0], params);

      if (imgParams.curTimeIdx == 0) {
        data->k0_time0[direction] = data->k0[direction];
      }
      /* check if the k0 vector found is reasonably close to the guess */
      vector deltak0;
      deltak0.x = data->k0[direction].x - data->k0guess[direction].x;
      deltak0.y = data->k0[direction].y - data->k0guess[direction].y;
      float dist = sqrt(deltak0.x * deltak0.x + deltak0.y * deltak0.y);
      if (dist > K0_WARNING_THRESH) {
        printf("WARNING: ");
      }
      printf("best fit for k0 is %f pixels from expected value.\n", dist);

      if (imgParams.ntimes > 1 && imgParams.curTimeIdx > 0
          && dist > 2*K0_WARNING_THRESH) {
        data->k0[direction] = data->k0_time0[direction];
        printf("k0 estimate of time point 0 is used instead\n");
        for (int order = 1; order < params->norders; ++order) {
          float corr_coeff;
          if (imgParams.nz0>1)
            corr_coeff = findrealspacemodamp(bands, &data->overlap0,
              &data->overlap1, imgParams.nx, imgParams.ny, imgParams.nz0,
              0, order, data->k0[direction], imgParams.dy, imgParams.dz,
              &data->otf, imgParams.wave[0], &data->amp[direction][order],
              &amp_inv, &amp_combo, 1, params);
          else
            corr_coeff = findrealspacemodamp(bands, &data->overlap0,
              &data->overlap1, imgParams.nx, imgParams.ny, imgParams.nz0,
              order-1, order, data->k0[direction], imgParams.dy, imgParams.dz,
              &data->otf, imgParams.wave[0], &data->amp[direction][order],
              &amp_inv, &amp_combo, 1, params);
          printf("modamp mag=%f, phase=%f\n, correlation coeff=%f\n\n",
                 cmag(data->amp[direction][order]),
                 atan2(data->amp[direction][order].y, data->amp[direction][order].x),
                 corr_coeff);
        }
      }
    } else {
      /* assume k0 vector known, so just fit for the modulation amplitude and phase */
      printf("known k0 for direction %d = (%f, %f) \n", direction, 
          data->k0[direction].x, data->k0[direction].y);
      for (int order = 1; order < params->norders; ++order) {
        float corr_coeff = findrealspacemodamp(bands, &data->overlap0,
            &data->overlap1, imgParams.nx, imgParams.ny, imgParams.nz0, 
            0, order, data->k0[direction], imgParams.dy, imgParams.dz,
            &data->otf, imgParams.wave[0], &data->amp[direction][order],
            &amp_inv, &amp_combo, 1, params);
        printf("modamp mag=%f, phase=%f\n",
            cmag(data->amp[direction][order]),
            atan2(data->amp[direction][order].y, data->amp[direction][order].x));
        printf("reverse modamp mag=%f, phase=%f\n", 1.0f / cmag(amp_inv),
            -atan2(amp_inv.y, amp_inv.x));
        printf("combined modamp mag=%f, phase=%f\n",
            cmag(amp_combo), atan2(amp_combo.y, amp_combo.x));
        printf("correlation coeff=%f\n\n", corr_coeff);
#ifndef __SIRECON_USE_TIFF__
        if (order == 1 && params->bSaveOverlaps) {// output the overlaps
          // output the overlaps
          CPUBuffer tmp0(data->overlap0.getSize());
          data->overlap0.set(&tmp0, 0, tmp0.getSize(), 0);
          CPUBuffer tmp1(data->overlap1.getSize());
          data->overlap1.set(&tmp1, 0, tmp1.getSize(), 0);
          float* ol0Ptr = (float*)tmp0.getPtr();
          float* ol1Ptr = (float*)tmp1.getPtr();
          for (int z = 0; z < imgParams.nz0; ++z) {
            IMWrSec(overlaps_stream_no, ol0Ptr + z * imgParams.nx * imgParams.ny);
            IMWrSec(overlaps_stream_no, ol1Ptr + z * imgParams.nx * imgParams.ny);
          }
        }
#endif
      }
    }     /* if(searchforvector) ... else ... */

    if (imgParams.nz == 1) {
      /* In 2D SIM, amp stores modamp's between each adjacent pair of
       * bands. We want to convert this to modamp w.r.t. order 0 */
      for (int order = 2; order < params->norders; ++order) {
        data->amp[direction][order] = cmul(data->amp[direction][order],
            data->amp[direction][order - 1]);
      }
    }

    if (params->forceamp[0] > 0.0) {
      /* force modamp's amplitude to be a value user provided (ideally
       * should be 1)  */
      for (int order = 1; order < params->norders; ++order) {
        float a = cmag(data->amp[direction][order]);
        if (a < params->forceamp[order-1]) {
          float ampfact = params->forceamp[order-1] / a;
          data->amp[direction][order].x *= ampfact;
          data->amp[direction][order].y *= ampfact;
          printf("modamp mag=%f, phase=%f  \n",
              cmag(data->amp[direction][order]),
              atan2(data->amp[direction][order].y, data->amp[direction][order].x));
        }
      }
    }

    /* In 2D NLSIM, we often don't trust the modamp fit between neighboring high-order components; */
    /* only if fitallphases is True do all fitted modamps get actually used; */
    /* otherwise, order 2 and above's global phase is inferred from order 1's phase.  */
    if (!params->bFitallphases) {
      float base_phase = get_phase(data->amp[direction][1]);
      cuFloatComplex expiphi;
      cuFloatComplex amplitude;
      float phi;
      for (int order = 2; order < params->norders; ++order) {
        amplitude.x = cmag(data->amp[direction][order]);
        amplitude.y = 0;
        phi = order * base_phase /*+ M_PI*(order-1)*/;
        /* sign flip for every other order happens only in saturation case? */
        expiphi.x = cos(phi);
        expiphi.y = sin(phi);
        data->amp[direction][order] = cmul(amplitude, expiphi);
      }
    }
  } /* end for (dir)  */   /* done finding the modulation vectors and phases for all directions */
}

void SIM_Reconstructor::loadImageData(int it, int iw, int zoffset)
{
#ifdef __SIRECON_USE_TIFF__
  // set up m_myParams, m_imgParams, and m_reconData based on the first input TIFF
  CImg<> rawtiff(m_all_matching_files[it].c_str());
  if (it == 0)
    ::setup(rawtiff, &m_myParams, &m_imgParams, &m_reconData);

  if (m_imgParams.nz0 < rawtiff.depth())
    rawtiff.crop(0,0,0,0,rawtiff.width()-1,rawtiff.height()-1,
                 m_imgParams.nz0*m_myParams.ndirs*m_myParams.nphases-1,0);
#endif

  for (int direction = 0; direction < m_myParams.ndirs; ++direction) {
    // Temporary Buffers for reading switch-off images
    PinnedCPUBuffer buffer(sizeof(float) * m_imgParams.nx * m_imgParams.ny);
    /*Pinned*/CPUBuffer offBuff(sizeof(float) * (m_imgParams.nx + 2) * m_imgParams.ny);

    std::vector<GPUBuffer>* rawImages = &(m_reconData.savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(m_reconData.savedBands[direction]);

    int z = 0;
    for (z = 0; z < m_imgParams.nz; ++z) {
      int zsec;
      /* First determine which section of the raw data to load */
      if (m_myParams.bFastSIM) {
        /* data organized into (nz, ndirs, nphases) */
        zsec = (z * m_myParams.ndirs * m_myParams.nphases +
                direction * m_myParams.nphases);
      } else { /* data organized into (ndirs, nz, nphases) */
        zsec = direction * m_imgParams.nz * m_myParams.nphases + z * m_myParams.nphases;
      }

      for (int phase = 0; phase < m_myParams.nphases; ++phase) {
#ifdef __SIRECON_USE_TIFF__
        load_and_flatfield(rawtiff, zsec, (float*)offBuff.getPtr(), m_myParams.constbkgd, 
                           m_imgParams.inscale);
#else
        if (m_myParams.bBgInExtHdr) {
          /* subtract the background value of each exposure stored in
           * extended header, indexed by the section number. */
          int extInts;
          float extFloats[3];
          IMRtExHdrZWT(istream_no, zsec, iw, it, &extInts, extFloats);
          m_reconData.backgroundExtra = extFloats[2];
        }
        load_and_flatfield(zsec, iw, it,
            (float*)offBuff.getPtr(), (float*)buffer.getPtr(),
            m_imgParams.nx, m_imgParams.ny,
            (float*)m_reconData.background.getPtr(), m_reconData.backgroundExtra,
            (float*)m_reconData.slope.getPtr(),
            m_imgParams.inscale, m_myParams.bUsecorr);
#endif
        assert(offBuff.hasNaNs() == false);
        // Transfer the data from offBuff to device buffer previously allocated (see m_reconData.savedBands)
        offBuff.set(&(rawImages->at(phase)),
            0, (m_imgParams.nx + 2) * m_imgParams.ny * sizeof(float),
            (z + zoffset) * (m_imgParams.nx + 2) * m_imgParams.ny * sizeof(float));

        ++zsec;
      } // end for (phase)
    } // end for (z)
  } // end for (direction)
}

void apodizationDriver(int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data)
{
  for (int direction = 0; direction < params->ndirs; ++direction) {
    /* data assumed taken with z changing fast, direction changing
     * slowly */
    std::vector<GPUBuffer>* rawImages = &(data->savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(data->savedBands[direction]);

    for (int z = 0; z < imgParams.nz; ++z) {
      for (int phase = 0; phase < params->nphases; ++phase) {
        if (params->napodize >= 0) {
          // Goes through here
          apodize(params->napodize, imgParams.nx, imgParams.ny,
              &(rawImages->at(phase)),
              (z + zoffset) * (imgParams.nx + 2) * imgParams.ny);
        } else if (params->napodize == -1) {
          cosapodize(imgParams.nx, imgParams.ny, &(*rawImages)[phase],
              (z + zoffset) * imgParams.nx * imgParams.ny);
        }
      } /* end for (phase), loading, flatfielding, and apodizing raw images */
    }
  } /* end for (dir)  */
}

void rescaleDriver(int it, int iw, int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data)
{
  for (int direction = 0; direction < params->ndirs; ++direction) {
    /* data assumed taken with z changing fast, direction changing
     * slowly */
    std::vector<GPUBuffer>* rawImages = &(data->savedBands[direction]);
    std::vector<GPUBuffer>* bands = &(data->savedBands[direction]);

    for (int z = 0; z < imgParams.nz; ++z) {
      if (params->do_rescale) {
        // Goes through here
        rescale(imgParams.nx, imgParams.ny, imgParams.nz, z, zoffset,
            direction, iw, it, params->nphases, rawImages, params->equalizez,
            params->equalizet, &data->sum_dir0_phase0[0]);
      }
    }
  }
}

void transformXYSlice(int zoffset, ReconParams* params,
    const ImageParams& imgParams, DriftParams* driftParams, ReconData* data)
{
  cufftHandle rfftplanGPU;
  int fftN[2] = {imgParams.ny, imgParams.nx};
  int inembed[2] = {imgParams.nx * imgParams.ny, imgParams.nx + 2};
  int istride = 1;
  int idist = (imgParams.nx + 2) * imgParams.ny;
  int onembed[2] = {imgParams.nx * imgParams.ny, imgParams.nx /2 +1};
  int ostride = 1;
  int odist = (imgParams.nx / 2 + 1) * imgParams.ny;
  cufftResult cuFFTErr = cufftPlanMany(&rfftplanGPU, 2, &fftN[0],
      inembed, istride, idist,
      onembed, istride, odist,
      CUFFT_R2C, imgParams.nz);
  if (cuFFTErr != CUFFT_SUCCESS) {
    std::cout << "Error code: " << cuFFTErr << std::endl;
    throw std::runtime_error("cufftPlanMany() failed.");
  }

  for (int direction = 0; direction < params->ndirs; ++direction) {
    for (int phase = 0; phase < params->nphases; ++phase) {
      cuFFTErr = cufftExecR2C(rfftplanGPU,
          (float*)data->savedBands[direction][phase].getPtr(),
          (cuFloatComplex*)data->savedBands[direction][phase].getPtr());
      if (cuFFTErr != CUFFT_SUCCESS) {
        std::cout << "Line:" << __LINE__ ;
        throw std::runtime_error("cufft failed.");
      }
    }
  }

  cufftDestroy(rfftplanGPU);
  if (cuFFTErr != CUFFT_SUCCESS) {
    std::cout << "Error code: " << cuFFTErr << std::endl;
    throw std::runtime_error("cufftDestroy() failed.");
  }
}

/***************************** makematrix ************************************/
/*     generates the matrix that is to be used to separate the raw indata    */
/*     into the different bands of sample information.                       */
/*  Two cases:
    1. If arrPhases is NULL, we're certain that phases are equally spaced within 2Pi range;
    2. Else, user supplies a list of phases actually used (e.g, when there is drift or 
    non-ideal SLM patterns are used), and noise variance amplification factors for each 
    order and direction is calculated based on the non-ideal separation matrix
    */
/*****************************************************************************/
void makematrix(int nphases, int norders, int dir, float *arrPhases,
    float *sepMatrix, float * noiseVarFactors)
{
  std::cout << "In makematrix." << std::endl;
  /*
     norders is not necessarily nphases/2+1, in cases such as nonlinear SIM
     where the number of output orders can be smaller.
     */
  int j, order;
  float phi, sqrt_nphases, sqrt_2;

  sqrt_nphases = sqrt((float)nphases);
  sqrt_2 = sqrt(2.);

  // norders could be less than (nphases+1)/2
  if (arrPhases == 0) { /* phase values equally spaced between 0 and 2Pi */
    phi = 2.0 * M_PI / nphases;
    for (j = 0; j < nphases; ++j) {
      sepMatrix[0 * nphases + j] = 1.0;
      for (order = 1; order < norders; ++order) {
        /* With this definition, bandplus = bandre + i bandim 
         * has coefficient exp() with unit normalization */
        sepMatrix[(2 * order - 1) * nphases + j] = cos(j * order * phi);
        sepMatrix[2 * order * nphases + j] = sin(j * order * phi);
      }
    }
  }
  else {/* User supplied phase values in arrPhases */
    int *ipvt, info, i, nCols;
    float *workspace, *forwardM, *A_t_A, unit=1., zero=0.;

    /* Here we differentiate between 2 options: */
    /* 1. use direct matrix inversion;          */
    /* 2. use least-square solution.            */

    nCols = 2 * norders - 1;

    if (nCols < nphases) {  // use least-square
      ipvt = (int*)malloc(nCols * sizeof(int));
      forwardM   = (float*)malloc(nphases*nCols*sizeof(float));
      A_t_A = (float*)malloc(nCols*nCols*sizeof(float));
      /* First construct the forward matrix*/
      for (i=0; i<nphases; i++) {
        forwardM[i] = 1.0 / nphases;
        for (order=1; order<norders; order++) {
          forwardM[i+(2*order-1)*nphases] = 2 * cos(order*arrPhases[i]) / nphases;
          forwardM[i+2*order*nphases]     = 2 * sin(order*arrPhases[i]) / nphases;
        }
      }

      /*       print_matrix("forward matrix:", nphases, nCols, forwardM, nphases); */

      /* multiply transposed forward matrix with forward matrix */
      sgemm_("T", "N", &nCols, &nCols, &nphases, &unit, forwardM,  &nphases,
          forwardM, &nphases,  &zero, A_t_A,  &nCols);
      /*       print_matrix("A transpose times A:", nCols, nCols, A_t_A, nCols); */

      /* Then invert the forward matrix to form the sep matrix*/
      sgetrf_(&nCols, &nCols, A_t_A, &nCols, ipvt, &info);
      if (info==0) { // successful factorization
        workspace = (float*)malloc(nCols*sizeof(float));
        sgetri_(&nCols, A_t_A, &nCols, ipvt, workspace, &nCols, &info);
        free(workspace);
      }
      else
        printf("WARNING: sgetri() returns non 0.\n");

      if (info!=0)
        printf("WARNING: sgetrf() returns non 0.\n");

      // multiply inverse A_t_A with transposed forward matrix
      sgemm_("N", "T", &nCols, &nphases,  &nCols, &unit, A_t_A,  &nCols,
          forwardM, &nphases,  &zero, sepMatrix, &nCols);

      /*       print_matrix("Separation Matrix:", nCols, nphases, sepMatrix, nCols); */
      /* transpose back to C-style matrix */
      matrix_transpose(sepMatrix, nphases, nCols);

      free(ipvt);
      free(forwardM);
      free(A_t_A);
    }
    else {  // use direct inversion
      ipvt = (int*)malloc(nphases*sizeof(int));
      /* First construct the forward matrix (in C-style storage convention ) */
      for (i=0; i<nphases; i++) {
        sepMatrix[i*nphases] = 1.0 / nphases;
        for (order=1; order<norders; order++) {
          sepMatrix[i*nphases+(2*order-1)] = 2 * cos(order*arrPhases[i]) / nphases;
          sepMatrix[i*nphases+2*order]     = 2 * sin(order*arrPhases[i]) / nphases;
        }
      }
      /* Then invert the forward matrix to form the sep matrix*/
      sgetrf_(&nphases, &nphases, sepMatrix, &nphases, ipvt, &info);
      if (info==0) { // successful factorization
        workspace = (float*)malloc(nphases*sizeof(float));
        sgetri_(&nphases, sepMatrix, &nphases, ipvt, workspace, &nphases, &info);
        free(workspace);
      }
      else
        printf("WARNING: sgetri() returns non 0.\n");

      if (info!=0)
        printf("WARNING: sgetrf() returns non 0.\n");
      free(ipvt);
    }

    /* Report noise factors */
    for (order=0;order<norders;order++)   noiseVarFactors[dir*norders+order] = 0.0;
    for (j=0;j<nphases;j++) {
      noiseVarFactors[dir*norders+0] += pow(sepMatrix[0*nphases+j], 2);
      for (order=1;order<norders;order++)
        noiseVarFactors[dir*norders+order] += pow(sepMatrix[(2*order-1)*nphases+j], 2) + pow(sepMatrix[2*order*nphases+j], 2);
    }

    printf(" Separation noise factors: ");
    for (order=0;order<norders;order++) {
      noiseVarFactors[dir*norders+order] /= nphases;
      printf(" Order %d: %.3f,",order, sqrt(noiseVarFactors[dir*norders+order]));
    }
    printf("\n");

  }

  printf("Separation matrix:\n");
  for (j=0;j<norders*2-1;j++) {
    int i;
    for (i=0; i<nphases; i++)
      printf("%9.5f ", sepMatrix[j*nphases+i]);
    printf("\n");
  }
  printf("\n");
}

void allocSepMatrixAndNoiseVarFactors(const ReconParams& params, ReconData* reconData)
{
  reconData->sepMatrix.resize(params.nphases * (params.norders * 2 - 1));
  reconData->noiseVarFactors.resize(params.ndirs * params.norders, 1.0f);
}

#ifdef __SIRECON_USE_TIFF__
void load_and_flatfield(CImg<> &cimg, int section_no, float *bufDestiny,
                        float background, float inscale)
{
  float *buffer = cimg.data(0, 0, section_no);
  int nx = cimg.width();
  int ny = cimg.height();
#pragma omp parallel for
  for (int l=0; l<ny; l++) {
    for (int k=0; k<nx; k++) {
      bufDestiny[l*(nx+2) + k] = (buffer[l*nx + k] - background) * inscale;
    }
    for(int k=nx;k<nx+2;k++)
      bufDestiny[l*(nx+2) + k] = 0.0;
  }
}
#else
void load_and_flatfield(int section_no, int wave_no, int time_no,
    float *bufDestiny, float *buffer, int nx, int ny, float *background,
    float backgroundExtra, float *slope, float inscale, int bUsecorr)
  /*
     Load the next 2D section from the MRC file identified by "istream_no".
     "bufDestiny" is where the current loaded section ends up being; it's assumed to have 2 extra columns for in-place FFT later
     "buffer" is a nx*ny sized array to hold temporarily the loaded data before it is flat-fielded and copied to "bufDestiny"
     */
{
  IMPosnZWT(istream_no, section_no, wave_no, time_no);
  IMRdSec(istream_no, buffer);

  if (bUsecorr) {
#pragma omp parallel for
    for (int l=0; l<ny; l++) {
      for (int k=0; k<nx; k++) {
        bufDestiny[l*(nx+2) + k] = ((buffer[l*nx + k]-background[l*nx+k]-backgroundExtra) * slope[l*nx+k]) * inscale;
      }
      for(int k=nx;k<nx+2;k++)
        bufDestiny[l*(nx+2) + k] = 0.0;
    }
  } else {
#pragma omp parallel for
    for (int l=0; l<ny; l++) {
      for (int k=0; k<nx; k++) {
        bufDestiny[l*(nx+2) + k] = (buffer[l*nx + k]-background[l*nx+k]-backgroundExtra) * inscale;
      }
      for(int k=nx;k<nx+2;k++)
        bufDestiny[l*(nx+2) + k] = 0.0;
    }
  }
}
#endif

void matrix_transpose(float* mat, int nRows, int nCols)
{
  int i, j;
  float* tmpmat = (float*)malloc(nRows*nCols*sizeof(float));

  for (i=0; i<nRows; i++)
    for (j=0; j<nCols; j++)
      tmpmat[j*nRows+i] = mat[i*nCols+j];
  memcpy(mat, tmpmat, nRows*nCols*sizeof(float));
  free(tmpmat);
}

int rdistcutoff(int iw, const ReconParams& params, const ImageParams& imgParams)
{
  float dkr = 1.0 / (imgParams.ny * imgParams.dy);
  int result = (int)floor((params.na * 2.0 / (
          imgParams.wave[iw] / 1000.0)) / dkr);
  return result;
}

int fitXYdrift(vector3d *drifts, float * timestamps, int nPoints, vector3d *fitted_drift, float *eval_timestamps, int nEvalPoints)
  /*
     fit a curve using the data points of (x,y) vectors in "drifts" versus "timestamps". Then evaluate drifts for
     the nEvalPoints points using the fitted curve; return the result in fitted_drifts.
     Note: if nPoints is greater than 3, the use cubic fit? otherwise use parabola fit?
     drifts[0] should always be (0, 0)?
     */
{
  int   i, j;
  int   nPolyTerms = 4; // highest order term is nPolyTerms-1; this could be passed in as an argument
  float *matrix_A, *vecRHS;
  float *workspace, wkopt;
  int   info, lwork, nRHS=2;  // 2 because of solving for both x and y
  int   bForceZero = 1;

  if (bForceZero) {  /* for the curve to go through time point 0 */
    int nPoints_minus1, nPolyTerms_minus1;
    /* subtract all time stamps by the first */
    for (i=1; i<nPoints; i++)
      timestamps[i] -= timestamps[0];

    /* Therefore, the solution to the constant is drifts[0].
       Now we have 2 coefficients to solve. */ 

    // construct the forward matrix; it's a transpose of the C-style matrix
    matrix_A = (float*)malloc((nPoints-1) * (nPolyTerms-1) * sizeof(float));
    // leading dimension is nPoints (column in Fortran, row in C)
    for (i=0; i<nPolyTerms-1; i++)
      for (j=0; j<nPoints-1; j++)
        matrix_A[i*(nPoints-1) + j] = pow(timestamps[j+1], i+1);

    // construct the right-hand-side vectors
    vecRHS = (float*)malloc((nPoints-1) * nRHS * sizeof(float));
    for (j=0; j<nPoints-1; j++) {
      vecRHS[j]           = drifts[j+1].x - drifts[0].x;
      vecRHS[j+nPoints-1] = drifts[j+1].y - drifts[0].y;
    }

    /* Because of calling Fortran subroutines, every parameter has to be passed as a pointer */
    nPoints_minus1 = nPoints - 1;
    nPolyTerms_minus1 = nPolyTerms - 1;
    /* Query and allocate the optimal workspace */
    lwork = -1;
    sgels_("No transpose", &nPoints_minus1, &nPolyTerms_minus1, &nRHS, matrix_A, &nPoints_minus1,
        vecRHS, &nPoints_minus1, &wkopt, &lwork, &info );
    lwork = wkopt;
    workspace = (float*)malloc( lwork* sizeof(float) );

    sgels_("No transpose", &nPoints_minus1, &nPolyTerms_minus1, &nRHS, matrix_A, &nPoints_minus1,
        vecRHS, &nPoints_minus1, workspace, &lwork, &info);

    if (info != 0)
      return 0;

    /* Now subtract timestamp[0] from all eval_timestamps */
    for (i=0; i<nEvalPoints; i++)
      eval_timestamps[i] -= timestamps[0];

    /* Evaluate at datapoints eval_timestamps */
    for (i=0; i<nEvalPoints; i++) {
      fitted_drift[i].x = drifts[0].x;
      fitted_drift[i].y = drifts[0].y;
      for (j=0; j<nPolyTerms-1; j++) {
        fitted_drift[i].x += vecRHS[j] * pow(eval_timestamps[i], j+1);
        fitted_drift[i].y += vecRHS[j+nPoints-1] * pow(eval_timestamps[i], j+1);
      }
    }
  }
  else {
    // construct the forward matrix; it's a transpose of the C-style matrix
    matrix_A = (float*)malloc(nPoints * nPolyTerms * sizeof(float));
    // leading dimension is nPoints (column in Fortran, row in C)
    for (i=0; i<nPolyTerms; i++)
      for (j=0; j<nPoints; j++)
        matrix_A[i*nPoints + j] = pow(timestamps[j], i);

    // construct the right-hand-side vectors
    vecRHS = (float*)malloc(nPoints * nRHS * sizeof(float));

    for (j=0; j<nPoints; j++) {
      vecRHS[j] = drifts[j].x;
      vecRHS[j+nPoints] = drifts[j].y;
    }

    /* Query and allocate the optimal workspace */
    lwork = -1;
    sgels_("No transpose", &nPoints, &nPolyTerms, &nRHS, matrix_A, &nPoints, vecRHS, &nPoints, &wkopt, &lwork, &info );
    lwork = wkopt;
    workspace = (float*)malloc( lwork* sizeof(float) );

    sgels_("No transpose", &nPoints, &nPolyTerms, &nRHS, matrix_A, &nPoints, vecRHS, &nPoints, workspace, &lwork, &info);

    if (info != 0)
      return 0;

    /* Evaluate at datapoints eval_timestamps */
    for (i=0; i<nEvalPoints; i++) {
      fitted_drift[i].x = 0;
      fitted_drift[i].y = 0;
      for (j=0; j<nPolyTerms; j++) {
        fitted_drift[i].x += vecRHS[j] * pow(eval_timestamps[i], j);
        fitted_drift[i].y += vecRHS[j+nPoints] * pow(eval_timestamps[i], j);
      }
    }
  }

  free(matrix_A);
  free(vecRHS);
  free(workspace);

  return 1;
}

void calcPhaseList(float * phaseList, vector3d *driftlist,
    float *phaseAbs, float k0angle, float linespacing, float dr,
    int nphases, int nz, int direction, int z)
  /*
     Calculate the actual pattern phase value for each exposure, taking into account drift estimation, acquisition-time
     phase values, and the k0 vector of this direction.
     phaseAbs -- acquisition-time absolute phase used, recorded in extended header
     */
{
  int p;
  float phistep, k0mag;
  vector3d kvec, drift;

  phistep = 2*M_PI / nphases;

  k0mag = 1/linespacing;
  kvec.x = k0mag * cos(k0angle);
  kvec.y = k0mag * sin(k0angle);

  /* The "standard" evenly-spaced phases + phase correction caused by drift correction - 
     acquisition-time phase correction */
  for (p=0; p<nphases; p++) {
    float phiDrift;
    drift = driftlist[direction*nz*nphases + z*nphases + p];
    phiDrift = (drift.x*kvec.x + drift.y*kvec.y) * dr * 2 * M_PI;

    if (phaseAbs && nz == 1)
      // here the absolute phases of the off images are used; not necessary but simply because phase starts at 0 for off images
      phaseList[p] = phaseAbs[direction*nphases + p]  + phiDrift;

    else
      phaseList[p] = p * phistep + phiDrift;

    printf("phaseAbs[%d]=%10.5f, phaseList[%d]=%10.5f, phase error=%8.3f, phase drift=%8.3f deg\n",
        p, phaseAbs[direction*nphases + p]*180/M_PI,
        p, phaseList[p]*180/M_PI, 
        (phaseList[p] - p *phistep)*180/M_PI,
        phiDrift*180/M_PI);
  }

}

float get_phase(cuFloatComplex b)
{
  return atan2(b.y, b.x);
}

float cmag(cuFloatComplex a)
{
  return (sqrt(a.x*a.x+a.y*a.y));
}

cuFloatComplex cmul(cuFloatComplex a, cuFloatComplex b)
{
  cuFloatComplex result;
  result.x = a.x * b.x - a.y * b.y;
  result.y = a.x * b.y + a.y * b.x;
  return result;
}

#ifndef __SIRECON_USE_TIFF__
void saveIntermediateDataForDebugging(const ReconParams& params)
{
  if (params.bSaveSeparated) {
    IMWrHdr(separated_stream_no,
        "separated bands of all directions", 1, 0, 1, 0);
    IMClose(separated_stream_no);
  }
  if (params.bSaveAlignedRaw) {
    IMWrHdr(aligned_stream_no,
        "drift-corrected raw images", 1, aligned_header.amin,
        aligned_header.amax, aligned_header.amean);
    IMClose(aligned_stream_no);
  }
  if (params.bSaveOverlaps) {
    IMWrHdr(overlaps_stream_no,
        "overlaps in real space", 1, 0, 1, 0);
    IMClose(overlaps_stream_no);
  }
  if (params.bSaveSeparated ||
      params.bSaveAlignedRaw ||
      params.bSaveOverlaps) {
    printf(
        "\nQuit processing because either bSaveSeparated, "
        "bSaveAlignedRaw, or bSaveOverlaps is TRUE\n");
    exit(0);
  }
}
#endif

#ifndef __SIRECON_USE_TIFF__
void saveCommandLineToHeader(int argc, char **argv, IW_MRC_HEADER &header)
{
  char titles[1000];
  titles[0] = '\0';
  for (int i = 3; i < argc; ++i) {
    strcat(titles, argv[i]);
    strcat(titles, " ");
  }
  IMAlLab(ostream_no, titles, strlen(titles) / 80 + 1);
  IMWrHdr(ostream_no, header.label, 1, header.amin, header.amax,
      header.amean);
}
#endif

void dumpBands(std::vector<GPUBuffer>* bands, int nx, int ny, int nz0)
{
  int n = 0;
  for (std::vector<GPUBuffer>::iterator i = bands->begin();
      i != bands->end(); ++i) {
    std::stringstream s;
#ifdef  __SIRECON_USE_TIFF__
    s << "band" << n << ".tif";
    TIFF * band_tiff = TIFFOpen(s.str().c_str(), "w");
    CPUBuffer buf(nx*ny*nz0*sizeof(float));
    i->set(&buf, 0, buf.getSize(), 0);
    float *ptr = (float*) buf.getPtr();
    for (int z=0; z<nz0; z++) {
      save_tiff(band_tiff, z, 0, 1, nx, ny, ptr, 0);
      ptr += nx * ny;
    }
    TIFFClose(band_tiff);
#endif
    std::stringstream ss;
    ss << "band" << n << ".dat";
    std::ofstream os(ss.str().c_str());
    i->dump(os, nx, 0, nx * ny * sizeof(float));
    ++n;
  }
}

#ifndef NDEBUG
#ifndef _WIN32
void deviceMemoryUsage()
{
  nvmlDeviceGetMemoryInfo(nvmldevice, &memoryStruct);
  printf("****** Used memory %lldM; free memory %lldM\n", memoryStruct.used/1024/1024, memoryStruct.free/1024/1024);
}
#endif
#endif

#include <ostream>
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " ")); 
  return os;
}

// IMPLEMENTATION OF CLASS SIM_RECONSTRUCTOR FOLLOWS:

SIM_Reconstructor::SIM_Reconstructor(int argc, char **argv)
{
  // Call this constructor from a command-line based program

#ifndef NDEBUG
#ifndef _WIN32
  nvmlInit();
  nvmlDeviceGetHandleByIndex(0, &nvmldevice);
#endif
#endif
  m_argc = argc;
  m_argv = argv;

  SetDefaultParams(&m_myParams);
  
  // define all the commandline and config file options
  setupProgramOptions();
  po::positional_options_description p;
  p.add("input-file", 1);
  p.add("output-file", 1);
  p.add("otf-file", 1);

  // parse the commandline
  store(po::command_line_parser(argc, argv).
        options(m_progopts).positional(p).run(), m_varsmap);

  if (m_varsmap.count("help")) {
    std::cout << m_progopts << "\n";
    exit(0);
  }

  notify(m_varsmap);

  if (m_config_file != "") {
    std::ifstream ifs(m_config_file.c_str());
    if (!ifs) {
      std::cout << "can not open config file: " << m_config_file << "\n";
      std::cout << "proceed without it\n";
    }
    else {
      // parse config file
      store(parse_config_file(ifs, m_progopts), m_varsmap);
      notify(m_varsmap);
    }
  }

   // fill in m_myParams fields that have not been set yet
  setParams();

  printf("nphases=%d, ndirs=%d\n", m_myParams.nphases, m_myParams.ndirs);
  
  // std::cout << m_myParams.bNoKz0 << std::endl;

#ifdef __SIRECON_USE_TIFF__
  /* Suppress "unknown field" warnings */
  TIFFSetWarningHandler(NULL);

  // gather all TIFF files with matching names under the same folder:
  m_all_matching_files = gatherMatchingFiles(std::string(m_myParams.ifiles), 
                                             std::string(m_myParams.ofiles));
  // In TIFF mode, m_myParams.ifiles refers to the name of the folder raw data resides in;
  // and m_myParams.ofiles refers to a pattern in all the raw data file names.
  m_imgParams.ntimes = m_all_matching_files.size();

#else
  /* Suppress IVE display of file headers */
  IMAlPrt(0);
#endif

  openFiles();
  // deviceMemoryUsage();

  m_zoffset = 0;
#ifndef __SIRECON_USE_TIFF__
  setup();
  ::setOutputHeader(m_myParams, m_imgParams, m_in_out_header);

  if (m_myParams.nzPadTo) {
    m_zoffset = (m_imgParams.nz0 - m_imgParams.nz) / 2;
  }

  bgAndSlope(m_myParams, m_imgParams, &m_reconData);
#endif
}

SIM_Reconstructor::~SIM_Reconstructor()
{
}

int SIM_Reconstructor::setupProgramOptions()
{

  // Declare a group of options that will be 
  // allowed both on command line and in a config file
  m_progopts.add_options()
    ("input-file", po::value<std::string>()->required(), "input file (or data folder in TIFF mode)")
    ("output-file", po::value<std::string>()->required(), "output file (or filename pattern in TIFF mode)")
    ("otf-file", po::value<std::string>()->required(), "OTF file")
    ("usecorr", po::value<std::string>(), "use the flat-field correction file provided")
    ("ndirs", po::value<int>(&m_myParams.ndirs)->default_value(3),
     "number of directions")
    ("nphases", po::value<int>(&m_myParams.nphases)->default_value(5),
     "number of phases per direction")
    ("nordersout", po::value<int>(&m_myParams.norders_output)->default_value(0),
     "number of output orders; must be <= norders")
    ("angle0", po::value<float>(&m_myParams.k0startangle)->default_value(1.648),
     "angle of the first direction in radians")
    ("ls", po::value<float>(&m_myParams.linespacing)->default_value(0.172),
     "line spacing of SIM pattern in microns")
    ("na", po::value<float>(&m_myParams.na)->default_value(1.2),
     "Detection numerical aperture")
    ("nimm", po::value<float>(&m_myParams.nimm)->default_value(1.33),
     "refractive index of immersion medium")
    ("zoomfact", po::value<float>(&m_myParams.zoomfact)->default_value(2.),
     "lateral zoom factor")
    ("explodefact", po::value<float>(&m_myParams.explodefact)->default_value(1.),
     "artificially exploding the reciprocal-space distance between orders by this factor")
    ("zzoom", po::value<int>(&m_myParams.z_zoom)->default_value(1),
     "axial zoom factor")
    ("nofilteroverlaps", po::value<int>(&m_myParams.bFilteroverlaps)->implicit_value(false),
     "do not filter the overlaping region between bands usually used in trouble shooting")
    ("background", po::value<float>(&m_myParams.constbkgd)->default_value(0.),
     "camera readout background")
    ("wiener", po::value<float>(&m_myParams.wiener)->default_value(0.01),
     "Wiener constant")
    ("forcemodamp", po::value< std::vector<float> >()->multitoken(),
     "modamps forced to these values")
    ("k0angles", po::value< std::string >(), //po::value< std::vector<float> >()->multitoken(),
     "user given pattern vector k0 angles for all directions")
    ("otfRA", po::value<int>(&m_myParams.bRadAvgOTF)->implicit_value(1),
     "using rotationally averaged OTF")
    ("fastSI", po::value<int>(&m_myParams.bFastSIM)->implicit_value(1),
     "SIM data is organized in Z->Angle->Phase order; default being Angle->Z->Phase")
    ("k0searchAll", po::value<int>(&m_myParams.bUseTime0k0)->implicit_value(0),
     "search for k0 at all time points")
    ("equalizez", po::value<int>(&m_myParams.equalizez)->implicit_value(true), 
     "bleach correcting for z")
    ("equalizet", po::value<int>(&m_myParams.equalizet)->implicit_value(true), 
     "bleach correcting for time")
    ("dampenOrder0", po::value<int>(&m_myParams.bDampenOrder0)->implicit_value(true),
     "dampen order-0 in final assembly")
    ("nosuppress", po::value<int>(&m_myParams.bSuppress_singularities)->implicit_value(false),
     "do not suppress DC singularity in final assembly (good idea for 2D/TIRF data)")
    ("nokz0", po::value<int>(&m_myParams.bNoKz0)->implicit_value(true),
     "do not use kz=0 plane of the 0th order in the final assembly")
    ("gammaApo", po::value<float>(&m_myParams.apoGamma)->default_value(1.f),
     "output apodization gamma; 1.0 means triangular apo")
    ("saveprefiltered", po::value<std::string>(),
     "save separated bands (half Fourier space) into a file and exit")
    ("savealignedraw", po::value<std::string>(),
     "save drift-fixed raw data (half Fourier space) into a file and exit")
    ("saveoverlaps", po::value<std::string>(),
     "save overlap0 and overlap1 (real-space complex data) into a file and exit")
    ("config,c", po::value<std::string>(&m_config_file)->default_value(""),
     "name of a file of a configuration.")
    ("2lenses", po::value<int>(&m_myParams.bTwolens)->implicit_value(1), "I5S data")
    ("help,h", "produce help message")
#ifdef __SIRECON_USE_TIFF__
    ("xyres", po::value<float>(&m_imgParams.dy)->default_value(0.1),
     "x-y pixel size (only used for TIFF files)")
    ("zres", po::value<float>(&m_imgParams.dz)->default_value(0.144),
     "z pixel size (only used for TIFF files)")
    ("wavelength", po::value<short>(&m_imgParams.wave[0])->default_value(530),
     "emission wavelength (only used for TIFF files)")
#endif
    ;

  return 0;
}

int SIM_Reconstructor::setParams()
{
  // "input-file", "output-file", and "otf-file" are now all required arguments.

  if (m_varsmap.count("input-file")) {
    strcpy(m_myParams.ifiles, m_varsmap["input-file"].as<std::string>().c_str());
    // m_myParams.ifilein = 1;
  }
  
  if (m_varsmap.count("output-file")) {
    strcpy(m_myParams.ofiles, m_varsmap["output-file"].as<std::string>().c_str());
    // m_myParams.ofilein = 1;
  }

  if (m_varsmap.count("otf-file")) {
    strcpy(m_myParams.otffiles, m_varsmap["otf-file"].as<std::string>().c_str());
    // m_myParams.otffilein = 1;
  }

  if (m_varsmap.count("usecorr")) {
    strcpy(m_myParams.corrfiles, m_varsmap["usecorr"].as<std::string>().c_str());
    m_myParams.bUsecorr = 1;
  }

  if (m_varsmap.count("forcemodamp")) {
    m_myParams.forceamp = m_varsmap["forcemodamp"].as< std::vector<float> >();
  }

  if (m_varsmap.count("wiener")) {
    if (m_myParams.wiener > 0)
      m_myParams.bUseEstimatedWiener = false;
    std::cout<< "wiener=" << m_myParams.wiener << std::endl;
  }

  if (m_varsmap.count("gammaApo")) {
    m_myParams.apodizeoutput = 2;
    std::cout << "gamma=" << m_myParams.apoGamma << std::endl;
  }

  if (m_varsmap.count("saveprefiltered")) {
    strcpy(m_myParams.fileSeparated, m_varsmap["saveprefiltered"].as<std::string>().c_str());
    m_myParams.bSaveSeparated = 1;
  }

  if (m_varsmap.count("savealignedraw")) {
    strcpy(m_myParams.fileRawAligned, m_varsmap["savealignedraw"].as<std::string>().c_str());
    m_myParams.bSaveAlignedRaw = 1;
  }

  if (m_varsmap.count("saveoverlaps")) {
    strcpy(m_myParams.fileOverlaps, m_varsmap["saveoverlaps"].as<std::string>().c_str());
    m_myParams.bSaveOverlaps = 1;
  }

  if (m_varsmap.count("k0angles")) {
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char> > tokens(m_varsmap["k0angles"].as< std::string >(), sep);
    for ( boost::tokenizer<boost::char_separator<char> >::iterator it = tokens.begin();
          it != tokens.end();
          ++it)
        m_myParams.k0angles.push_back(strtod(it->c_str(), NULL));
    std::cout << m_myParams.k0angles << std::endl;
  }

  return 0;
}

void SIM_Reconstructor::openFiles()
{
#ifdef __SIRECON_USE_TIFF__
  if (!m_all_matching_files.size()) // TIFF files are not opened till loadAndRescaleImage()
#else
  if (IMOpen(istream_no, m_myParams.ifiles, "ro"))
#endif
    throw std::runtime_error("Input file not found");

  /* Create output file */
  // In TIFF mode, output files are not created until writeResult() is called
#ifndef __SIRECON_USE_TIFF__
  if (IMOpen(ostream_no, m_myParams.ofiles, "new")) {
    std::cerr << "File " << m_myParams.ofiles << " can not be created.\n";
    throw std::runtime_error("File not found");
  }
#endif

#ifdef __SIRECON_USE_TIFF__
  if (!(otf_tiff = TIFFOpen(m_myParams.otffiles, "r")))
  // m_otf_tiff.assign(m_myParams.otffiles); // will throw CImgIOException if file cannot be opened
#else
  if (IMOpen(otfstream_no, m_myParams.otffiles, "ro"))
#endif
    throw std::runtime_error("OTF file not found");
}



void SIM_Reconstructor::processOneVolume()
{
  // process one SIM volume (i.e., for the time point timeIdx)

  int zoffset = 0;
  if (m_myParams.nzPadTo) {
    zoffset = (m_imgParams.nz0 - m_imgParams.nz) / 2;
  }

  m_reconData.bigbuffer.resize(0);
  m_reconData.outbuffer.resize(0);

  m_reconData.overlap0.resize(m_imgParams.nx * m_imgParams.ny * m_imgParams.nz *
      sizeof(cuFloatComplex));
  m_reconData.overlap0.setToZero();
  m_reconData.overlap1.resize(m_imgParams.nx * m_imgParams.ny * m_imgParams.nz *
      sizeof(cuFloatComplex));
  m_reconData.overlap1.setToZero();

  findModulationVectorsAndPhasesForAllDirections(zoffset,
      &m_myParams, m_imgParams, &m_driftParams, &m_reconData);

  m_reconData.overlap0.resize(0);
  m_reconData.overlap1.resize(0);

  // deviceMemoryUsage();

#ifndef __SIRECON_USE_TIFF__
  saveIntermediateDataForDebugging(m_myParams);
#endif
  m_reconData.bigbuffer.resize((m_myParams.zoomfact * m_imgParams.nx) *
      (m_myParams.zoomfact * m_imgParams.ny) * (m_myParams.z_zoom * m_imgParams.nz0) *
      sizeof(cuFloatComplex));
  m_reconData.bigbuffer.setToZero();
  m_reconData.outbuffer.resize((m_myParams.zoomfact * m_imgParams.nx) *
      (m_myParams.zoomfact * m_imgParams.ny) * (m_myParams.z_zoom * m_imgParams.nz0) *
      sizeof(float));
  m_reconData.outbuffer.setToZero();

  // deviceMemoryUsage();

  for (int direction = 0; direction < m_myParams.ndirs; ++direction) {

    filterbands(direction, &m_reconData.savedBands[direction],
        m_reconData.k0, m_myParams.ndirs, m_myParams.norders,
        m_reconData.otf, m_imgParams.dy, m_imgParams.dz,
        m_reconData.amp, m_reconData.noiseVarFactors,
        m_imgParams.nx, m_imgParams.ny, m_imgParams.nz0, m_imgParams.wave[0],
        &m_myParams);
    assemblerealspacebands(direction, &m_reconData.outbuffer,
        &m_reconData.bigbuffer, &m_reconData.savedBands[direction],
        m_myParams.ndirs, m_myParams.norders, m_reconData.k0,
        m_imgParams.nx, m_imgParams.ny, m_imgParams.nz0,
        m_myParams.zoomfact, m_myParams.z_zoom, m_myParams.explodefact);
  }
}

#ifdef __SIRECON_USE_TIFF__
void SIM_Reconstructor::setup(CImg<> &inTIFF)
{
  m_imgParams.nx = inTIFF.width();
  m_imgParams.ny = inTIFF.height();
  m_imgParams.nz = inTIFF.depth();
  m_imgParams.nwaves = inTIFF.spectrum();
  m_imgParams.nz /= m_imgParams.nwaves * params->nphases * params->ndirs;
  //  dy, dz, wavelength[0] from command line input or config file
  if (params->nzPadTo) {
    m_imgParams.nz0 = params->nzPadTo;
  } else {
    m_imgParams.nz0 = m_imgParams.nz;
  }

  printf("nx=%d, ny=%d, nz=%d, nz0 = %d, nwaves=%d, ntimes=%d\n",
      m_imgParams.nx, m_imgParams.ny, m_imgParams.nz, m_imgParams.nz0,
      m_imgParams.nwaves, m_imgParams.ntimes);

  ::setup_part2(&m_myParams, &m_imgParams, &m_reconData);

}
#else
void SIM_Reconstructor::setup()
{
  ::loadHeader(m_myParams, &m_imgParams, m_in_out_header);
  ::setup_part2(&m_myParams, &m_imgParams, &m_reconData);
}
#endif

void SIM_Reconstructor::loadAndRescaleImage(int timeIdx, int waveIdx)
{
  loadImageData(timeIdx, waveIdx, m_zoffset);

  ::rescaleDriver(timeIdx, waveIdx, m_zoffset, &m_myParams, m_imgParams, 
                    &m_driftParams, &m_reconData);
}

void SIM_Reconstructor::writeResult(int it, int iw)
{
  CPUBuffer outbufferHost(
      (m_myParams.zoomfact * m_imgParams.nx) *
      (m_myParams.zoomfact * m_imgParams.ny) *
      (m_myParams.z_zoom * m_imgParams.nz0) *
      sizeof(float));
  // if (it==0)
  //   computeAminAmax(&reconData.outbuffer, m_myParams.zoomfact * m_imgParams.nx,
  //                   m_myParams.zoomfact * m_imgParams.ny, m_myParams.z_zoom * m_imgParams.nz,
  //                   &minval, &maxval);
  m_reconData.outbuffer.set(&outbufferHost, 0, outbufferHost.getSize(), 0);

#ifndef __clang__
  float t1 = omp_get_wtime();
#endif

#ifdef __SIRECON_USE_TIFF__
  CImg<> outCimg((float*) outbufferHost.getPtr(),
    m_myParams.zoomfact * m_imgParams.nx,
    m_myParams.zoomfact * m_imgParams.ny,
    m_myParams.z_zoom * m_imgParams.nz0, true);  // "true" means outCimg does not allocate host memory

  outCimg.save(makeOutputFilePath(m_all_matching_files[it], std::string("_proc")).c_str());

#else

  int zoffset = 0;
  if (m_myParams.nzPadTo) {
    zoffset = (m_imgParams.nz0 - m_imgParams.nz) / 2;
  }
  float* ptr = ((float*)outbufferHost.getPtr()) +
    (int)(zoffset * m_myParams.z_zoom * m_imgParams.nx * m_imgParams.ny *
        m_myParams.zoomfact * m_myParams.zoomfact);
  for (int i = 0; i < m_imgParams.nz * m_myParams.z_zoom; ++i) {
    IMWrSec(ostream_no, ptr);
    if (it == 0) {
      for (int j = 0;
          j < (int)(m_imgParams.nx * m_imgParams.ny *
            m_myParams.zoomfact * m_myParams.zoomfact);
          ++j) {
        if (ptr[j] > maxval) {
          maxval = ptr[j];
        } else if (ptr[j] < minval) {
          minval = ptr[j];
        }
      }
    }
    ptr += (int)(m_myParams.zoomfact * m_imgParams.nx *
        m_myParams.zoomfact * m_imgParams.ny);
  }
  if (it == 0 && iw == 0) {
    m_in_out_header.amin = minval;
    m_in_out_header.amax = maxval;
  }
#endif

#ifndef __clang__
  float t2 = omp_get_wtime();
  printf("amin, amax took: %f s\n", t2 - t1);
#endif

  printf("Time point %d, wave %d done\n", it, iw);
}

void SIM_Reconstructor::closeFiles()
{
#ifndef __SIRECON_USE_TIFF__
  ::IMClose(istream_no);
  ::saveCommandLineToHeader(m_argc, m_argv, m_in_out_header);
  ::IMClose(ostream_no);
#endif
}
