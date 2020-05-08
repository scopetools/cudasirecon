#ifndef GPU_FUNCTIONS_IMPL_H
#define GPU_FUNCTIONS_IMPL_H

#include "gpuFunctions.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../cutilSafeCall.h"

#ifndef RED_BLOCK_SIZE_X
#define RED_BLOCK_SIZE_X 32
#endif
#ifndef RED_BLOCK_SIZE_Y
#define RED_BLOCK_SIZE_Y 32
#endif


/** Data from pParams going into constant memory */
__constant__ int const_pParams_bSuppress_singularities;
__constant__ int const_pParams_suppression_radius;
__constant__ int const_pParams_bDampenOrder0;
__constant__ int const_pParams_bNoKz0;
__constant__ int const_pParams_bFilteroverlaps;
__constant__ int const_pParams_apodizeoutput;
__constant__ float const_pParams_apoGamma;
__constant__ int const_pParams_bBessel;
__constant__ int const_pParams_bRadAvgOTF;
__constant__ int const_pParams_nzotf;
__constant__ float const_pParams_dkrotf;
__constant__ float const_wiener;

/** These data are not modified in the kernels and can go in constant
 * memory */
__constant__ int const_zdistcutoff[3];
__constant__ float const_ampmag2[3];
__constant__ float const_ampmag2_alldirs[9];
__constant__ cuFloatComplex const_conjamp[3];
__constant__ float const_noiseVarFactors[9];
__constant__ float2 const_k0[3];
__constant__ cuFloatComplex * const_otfPtrs[3];

#define MAX_ORDERS 32
#define MAX_PHASES 32
__constant__ float* const_outputPtrs[MAX_ORDERS * 2 -1];
__constant__ float* const_imgPtrs[MAX_PHASES];
__constant__ float const_sepMatrix[(MAX_ORDERS * 2 - 1) * MAX_PHASES];

__global__ void image_arithmetic_kernel(float* a, const float* b,
    int len, float alpha, float beta);
__global__ void apodize_x_kernel(int napodize, int nx, int ny,
    float* image);
__global__ void apodize_y_kernel(int napodize, int nx, int ny,
    float* image);
__global__ void cosapodize_kernel(int nx, int ny, float* image);
__global__ void rescale_kernel(float* img, int nx, int ny,
    float scaleFactor);
__global__ void sum_reduction_kernel(float* img, int nx, int ny,
    float* partialReduction);

__host__ void makeoverlaps(std::vector<GPUBuffer>* bands,
    GPUBuffer* overlap0, GPUBuffer* overlap1, int nx, int ny, int nz,
    int order1, int order2,
    float k0x, float k0y, float dy, float dz,
    std::vector<GPUBuffer>* OTF, short wave, ReconParams* pParams);

__global__ void makeOverlaps0Kernel(int nx, int ny, int nz,
    int order1, int order2, float kx, float ky,
    float rdistcutoff, float otfcutoff, float zdistcutoff,
	float order0_2_factor, float dkx, float dky, float kzscale,
    cuFloatComplex *band1im, cuFloatComplex *band1re,
    cuFloatComplex *overlap0);
__global__ void makeOverlaps1Kernel(int nx, int ny, int nz,
    int order1, int order2, float kx, float ky,
    float rdistcutoff, float otfcutoff, float zdistcutoff,
	float order0_2_factor, float dkx, float dky, float kzscale,
    cuFloatComplex *band2im, cuFloatComplex *band2re,
    cuFloatComplex *overlap1);

__host__ void aTimesConjB(GPUBuffer* overlap0, GPUBuffer* overlap1,
    int nx, int ny, int nz, GPUBuffer* crosscorr_c);
__global__ void aTimesConjBKernel(cuFloatComplex* overlap0,
    cuFloatComplex* overlap1, int nx, int ny, int nz,
    cuFloatComplex* crosscorr_c);

__host__ void computeIntensities(GPUBuffer* amplitudes, int nx, int ny,
    GPUBuffer* intensities);
__global__ void computeIntensitiesKernel(cuFloatComplex* amplitudes,
    int nx, int ny, float* intensities);

__host__ void findpeak(float array[], int sizex, int sizey, vector *peak);
__host__ float fitparabola( float a1, float a2, float a3);

__host__ void computeIntensities(GPUBuffer* amplitudes, int nx, int ny,
    GPUBuffer* intensities);

__host__ float getmodamp(float kangle, float klength,
    std::vector<GPUBuffer>* bands, GPUBuffer* overlap0, GPUBuffer* overlap1,
    int nx, int ny,int nz, int order1, int order2, float dy, float dz,
    std::vector<GPUBuffer>* otf, short wave, cuFloatComplex* modamp,
    int redoarrays, ReconParams *pParams, int bShowDetail); 

__host__ float findrealspacemodamp(std::vector<GPUBuffer>* bands,
    GPUBuffer* overlap0, GPUBuffer* overlap1, int nx, int ny, int nz,
    int order1, int order2, vector k0, float dy, float dz,
    std::vector<GPUBuffer>* OTF, short wave, cuFloatComplex* modamp1,
    cuFloatComplex* modamp2, cuFloatComplex* modamp3, int redoarrays,
    ReconParams *pParams);

__global__ void reductionKernel(
    int nx, int ny, int nz,
    float kx, float ky, float dxy,
    const cuFloatComplex *overlap0, const cuFloatComplex *overlap1,
    cuFloatComplex *XStarY, float *sumXMag, float *sumYMag);

__host__ float fitxyparabola( float x1, float y1, float x2, float y2,
    float x3, float y3);

__device__ cuFloatComplex dev_otfinterpolate(cuFloatComplex * otf, float
    kx, float ky, int kz, float kzscale);

__device__ float dev_suppress(float x);
__device__ float dev_mag2(cuFloatComplex x);
__device__ float dev_order0damping(float radius, float zindex, float rlimit, int zlimit);
__global__ void move_kernel(cuFloatComplex *inarray1, cuFloatComplex *inarray2, int order, 
			    cuFloatComplex *outarray, int nx, int ny, int nz, float
          zoomfact, int z_zoom);
__global__ void write_outbuffer_kernel1(cuFloatComplex * bigbuffer,
                                        float * outbuffer, int);
__global__ void write_outbuffer_kernel2(float * coslookup, float * sinlookup, 
                                        cuFloatComplex * bigbuffer, float * outbuffer, int);

__global__ void cos_sin_kernel(float k0x, float k0y, float dxy, float fact,
							   float * coslookup, float * sinlookup, int);

__global__ void filterbands_kernel1(int dir, int ndirs, int order, int norders, int nx, 
    int ny, int nz, float rdistcutoff, float zapocutoff, float apocutoff, 
	float krscale, float kzscale,
    cuFloatComplex * dev_bandptr, cuFloatComplex * dev_bandptr2, bool bSecondEntry);

// __global__ void filterbands_kernel2(int dir, int ndirs, int order, int norders, int nx, 
//     int ny, int nz, float rdistcutoff, float zapocutoff, float apocutoff, 
// 	float krscale, float kzscale, /*cuFloatComplex *dev_scale,*/ cuFloatComplex * dev_tempbandplus,
//     cuFloatComplex * dev_bandptr, cuFloatComplex * dev_bandptr2);

__global__ void filterbands_kernel3(int order, int nx, int ny, int nz,
				    cuFloatComplex * dev_bandptr, cuFloatComplex *
            dev_bandptr2);

__global__ void filterbands_kernel4(int order, int nx, int ny, int nz,
    cuFloatComplex * dev_tempbandplus, cuFloatComplex * dev_bandptr,
    cuFloatComplex * dev_bandptr2);

__global__ void separate_kernel(int norders, int nphases, int nx, int ny, int nz);
__global__ void computeAminAmax_kernel(const float* data, int numElems,
    float* maxPartialResult, float* minPartialResult);

__global__ void summation_kernel(float * img, double * intRes, int n);
__global__ void sumAboveThresh_kernel(float * img, double * intRes, unsigned * counter, float thresh, int n);
__global__ void scale_kernel(float * img, double factor, int n);


template<typename T>
T cpuReduce(const T* vec, int n) {
  T red = 0;
  for (int i = 0; i < n; ++i) {
    red += vec[i];
  }
  return red;
}

cuFloatComplex cpuReduce(const cuFloatComplex* vec, int n) {
  cuFloatComplex red;
  red.x = 0;
  red.y = 0;
  for (int i = 0; i < n; ++i) {
    red.x += vec[i].x;
    red.y += vec[i].y;
  }
  return red;
}
#endif
