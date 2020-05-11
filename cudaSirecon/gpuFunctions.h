#ifndef GPU_FUNCTIONS_H
#define GPU_FUNCTIONS_H

#include "GPUBuffer.h"
#include "cudaSireconImpl.h"

#include <vector>

struct vector;
struct vector3d;
struct ReconParams;

/*
  returns alpha*a + beta*b in a
  offset is an offset into a
  len is the number of elements in the vectors a and b
  imSrc can be a NULL pointer, in which case only fDest*imDest will be returned
*/
void image_arithmetic(GPUBuffer* a, const GPUBuffer& b, int offset,
    int len, float alpha, float beta);
void image_arithmetic(GPUBuffer* a, const GPUBuffer& b,
    int offsetA, int offsetB,
    int len, float alpha, float beta);

void apodize(int napodize, int nx,int ny, GPUBuffer* image, int offset);
void cosapodize(int nx,int ny, GPUBuffer* image, int offset);
void rescale(int nx, int ny, int nz, int z, int zoffset, int direction,
    int wave, int t, int phases, std::vector<GPUBuffer>* images, int equalizez,
    int equalizet, double* sum_dir0_phase0);

int calcRefImage(const std::vector<GPUBuffer>& rawImages,
    GPUBuffer* refImage, const std::vector<GPUBuffer>& offImages,
    int nOffImages, int nx, int ny, int nphases, int type_of_refImage);

void fixdrift_2D(std::vector<GPUBuffer>* CrawImages,
    vector3d *driftlist, int nphases, int nx, int ny, int nz, int dir,
    int z);

void separate(int nx, int ny, int z, int direction, int nphases, int
    norders, std::vector<GPUBuffer>*rawImages, float *sepMatrix);

void makemodeldata(int nx, int ny, int nz, std::vector<GPUBuffer>* bands,
    int norders, vector k0, float dy, float dz,
    std::vector<GPUBuffer>* OTF, short wave, ReconParams *pParams);

void fixdrift_bt_dirs(std::vector<GPUBuffer>* bands, int norders, 
    vector3d drift, int nx,int ny, int nz);

void findk0(std::vector<GPUBuffer>* bands, GPUBuffer* overlap0,
    GPUBuffer* overlap1, int nx, int ny, int nz, int norders, vector *k0,
    float dy, float dz, std::vector<GPUBuffer>* OTF, short wave,
    ReconParams * pParams);

void fitk0andmodamps(std::vector<GPUBuffer>* bands, GPUBuffer* overlap0,
    GPUBuffer* overlap1, int nx, int ny, int nz, int norders,
    vector *k0, float dxy, float dz, std::vector<GPUBuffer>* otf, short wave,
    cuFloatComplex* amps, ReconParams * pParams);

float findrealspacemodamp(std::vector<GPUBuffer>* bands,
    GPUBuffer* overlap0, GPUBuffer* overlap1, int nx, int ny, int nz,
    int order1, int order2, vector k0, float dy, float dz,
    std::vector<GPUBuffer>* OTF, short wave, cuFloatComplex* modamp1,
    cuFloatComplex* modamp2, cuFloatComplex* modamp3, int redoarrays,
    ReconParams *pParams);

void filterbands(int dir, std::vector<GPUBuffer>* bands,
    const std::vector<vector>& k0, int ndirs, int norders,
    std::vector<GPUBuffer>& otf, float dy, float dz, 
    const std::vector<std::vector<cuFloatComplex> >& amp,
    const std::vector<float>& noiseVarFactors, int nx, int ny, int nz,
    short wave, ReconParams* params);

void assemblerealspacebands(int dir, GPUBuffer* outbuffer, GPUBuffer* bigbuffer,
    std::vector<GPUBuffer>* bands, int ndirs, int norders,
    const std::vector<vector>& k0, int nx, int ny, int nz,
    float dxy, float zoomfact,
    int z_zoom, float expfact);

void computeAminAmax(const GPUBuffer* data, int nx, int ny, int nz,
    float* min, float* max);

#endif
