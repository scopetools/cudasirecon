#include "gpuFunctionsImpl_hh.cu"

__host__ void image_arithmetic(GPUBuffer* a, const GPUBuffer& b, int offset,
    int len, float alpha, float beta)
{
  float* aPtr = (float*)a->getPtr();
  aPtr += offset;
  const float* bPtr = (const float*)b.getPtr();

  int blockSize = 128;
  int numBlocks = (int)(ceil((float)len / blockSize));
  image_arithmetic_kernel<<<numBlocks, blockSize>>>(aPtr, bPtr,
      len, alpha, beta);
}

__host__ void image_arithmetic(GPUBuffer* a, const GPUBuffer& b,
    int offsetA, int offsetB, int len, float alpha, float beta)
{
  float* aPtr = (float*)a->getPtr();
  aPtr += offsetA;
  const float* bPtr = (const float*)b.getPtr();
  bPtr += offsetB;

  int blockSize = 128;
  int numBlocks = (int)(ceil((float)len / blockSize));
  image_arithmetic_kernel<<<numBlocks, blockSize>>>(aPtr, bPtr,
      len, alpha, beta);
}

__global__ void image_arithmetic_kernel(float* a, const float* b,
    int len, float alpha, float beta)
{
  int tid = blockDim.x * blockIdx.x + threadIdx.x;
  if (tid < len) {
    a[tid] = alpha * a[tid] + beta * b[tid];
  }
}

__host__ void apodize(int napodize, int nx,int ny, GPUBuffer* image,
    int offset)
{
  int blockSize = 64;
  int numBlocks = (int)(ceil((float)nx / blockSize));
  apodize_x_kernel<<<numBlocks, blockSize>>>(napodize, nx, ny,
      ((float*)image->getPtr()) + offset);

  numBlocks = (int)(ceil((float)ny / blockSize));
  apodize_y_kernel<<<numBlocks, blockSize>>>(napodize, nx, ny,
      ((float*)image->getPtr()) + offset);
}

__global__ void apodize_x_kernel(int napodize, int nx, int ny,
    float* image)
{
  int k = blockDim.x * blockIdx.x + threadIdx.x;
  if (k < nx) {
    float diff = (image[(ny - 1) * (nx + 2) + k] - image[k]) / 2.0;
    for (int l = 0; l < napodize; ++l) {
      float fact = 1.0 - sin((((float)l + 0.5) / (float)napodize) *
          M_PI * 0.5);
      image[l * (nx + 2) + k] = image[l * (nx + 2) + k] + diff * fact;
      image[(ny - 1 - l) * (nx + 2) + k] = image[(ny - 1 - l) * (nx + 2) + k] -
        diff * fact;
    }
  }
}

__global__ void apodize_y_kernel(int napodize, int nx, int ny,
    float* image)
{
  int l = blockDim.x * blockIdx.x + threadIdx.x;
  if (l < ny) {
    float diff = (image[l * (nx + 2) + nx - 1] - image[l * (nx + 2)]) / 2.0;
    for (int k = 0; k < napodize; ++k) {
      float fact = 1.0 - sin(((k + 0.5) / (float)napodize) * M_PI *
          0.5);
      image[l * (nx + 2) + k] = image[l * (nx + 2) + k] + diff * fact;
      image[l * (nx + 2) + (nx - 1 - k)] =
        image[l * (nx + 2) + (nx - 1 - k)] - diff * fact;
    }
  }
}

__host__ void cosapodize(int nx,int ny, GPUBuffer* image, int offset)
{
  dim3 blockSize;
  blockSize.x = 16;
  blockSize.y = 16;
  blockSize.z = 1;
  dim3 numBlocks;
  numBlocks.x = (int)(ceil((float)nx / blockSize.x));
  numBlocks.y = (int)(ceil((float)ny / blockSize.y));
  numBlocks.z = 1;
  cosapodize_kernel<<<numBlocks, blockSize>>>(nx, ny,
      ((float*)image->getPtr()) + offset);
}

__global__ void cosapodize_kernel(int nx, int ny, float* image)
{
  int k = blockDim.x * blockIdx.x + threadIdx.x;
  int l = blockDim.y * blockIdx.y + threadIdx.y;
  if (k<nx && l<ny) {
    float xfact = sin(M_PI * ((float)k + 0.5) / nx);
    float yfact = sin(M_PI * ((float)l + 0.5) / ny);
    image[l * (nx + 2) + k] *= xfact * yfact;
  }
}

__host__ void rescale(int nx, int ny, int nz, int z, int zoffset, int direction,
    int wave, int t, int nphases, std::vector<GPUBuffer>* images, int equalizez,
    int equalizet, double* sum_dir0_phase0)
{
  std::vector<float> sum(nphases);

  dim3 blockSize;
  blockSize.x = RED_BLOCK_SIZE_X;
  blockSize.y = RED_BLOCK_SIZE_Y;
  blockSize.z = 1;
  dim3 numBlocks;
  numBlocks.x = (int)(ceil((float)nx / blockSize.x));
  numBlocks.y = (int)(ceil((float)ny / blockSize.y));
  numBlocks.z = 1;

  GPUBuffer sumTmpDev(numBlocks.x * numBlocks.y * sizeof(float), 0);
  CPUBuffer sumTmpHost(numBlocks.x * numBlocks.y * sizeof(float));
  for (int phase = 0; phase < nphases; ++phase) {
    sum_reduction_kernel<<<numBlocks, blockSize>>>(
         ((float*)(images->at(phase).getPtr())) +
         (z + zoffset) * (nx + 2) * ny, nx, ny,
        (float*)sumTmpDev.getPtr());
    sumTmpDev.set(&sumTmpHost, 0, sumTmpDev.getSize(), 0);
    sum[phase] = 0.0f;
    for (int i = 0; i < numBlocks.x * numBlocks.y; ++i) {
      sum[phase] += ((float*)sumTmpHost.getPtr())[i];
    }
  }

  if (direction == 0 && !(equalizet && t != 0)) {
    sum_dir0_phase0[wave * nz + z] = (double)sum[0];
  }
  float ref;
  if (equalizez) {
    ref = sum_dir0_phase0[wave * nz + 0];
  } else  {
    ref = sum_dir0_phase0[wave * nz + z];
  }

  for (int phase = 0; phase < nphases; ++phase) {
    float ratio = ref / sum[phase];
    rescale_kernel<<<numBlocks, blockSize>>>(
         (float*)((*images)[phase].getPtr()) +
         (z + zoffset) * (nx + 2) * ny, nx, ny,
        ratio);
  }
}

__global__ void sum_reduction_kernel(float* img, int nx, int ny,
    float* partialReduction)
{
  int k = blockDim.x * blockIdx.x + threadIdx.x;
  int l = blockDim.y * blockIdx.y + threadIdx.y;
  __shared__ float locRedBuffer[RED_BLOCK_SIZE_X * RED_BLOCK_SIZE_Y];
  if (k < nx && l < ny) {
    locRedBuffer[threadIdx.x + blockDim.x * threadIdx.y] =
      img[k + l * (nx+2)];
  } else {
    locRedBuffer[threadIdx.x + blockDim.x * threadIdx.y] = 0.0f;
  }
  int ltid = threadIdx.y * RED_BLOCK_SIZE_X + threadIdx.x;
  __syncthreads();
  for (int s = RED_BLOCK_SIZE_X * RED_BLOCK_SIZE_Y / 2; s > 0; s >>= 1) {
    if (ltid < s) {
      locRedBuffer[ltid] += locRedBuffer[ltid + s];
    }
    __syncthreads();
  }
  if (ltid == 0) {
    int blockIndex = blockIdx.y * gridDim.x + blockIdx.x;
    partialReduction[blockIndex] = locRedBuffer[0];
  }
}

__global__ void rescale_kernel(float* img, int nx, int ny,
    float scaleFactor)
{
  int k = blockDim.x * blockIdx.x + threadIdx.x;
  int l = blockDim.y * blockIdx.y + threadIdx.y;
  if (k < nx && l < ny) {
    img[l * (nx + 2) + k] *= scaleFactor;
  }
}

__host__ float estimate_Wiener(const std::vector<GPUBuffer>& rawImages, int nx,
          int ny, int z, int nphases, int rdistcutoff)
{
  printf("In estimate_Wiener.\n");
  fflush(stdout);
  return 0.0f;
}

__host__ void fixdrift_2D(std::vector<GPUBuffer>* CrawImages,
    vector3d *driftlist, int nphases, int nx, int ny, int nz, int dir,
    int z)
{
  printf("In fixdrift_2D.\n");
  fflush(stdout);
}

__host__ int calcRefImage(const std::vector<GPUBuffer>& rawImages,
    GPUBuffer* refImage, const std::vector<GPUBuffer>& offImages,
    int nOffImages, int nx, int ny, int nphases, int type_of_refImage)
{
  printf("In calcRefImage.\n");
  fflush(stdout);
  return 0;
}

__host__ void determinedrift_2D(const std::vector<GPUBuffer>& rawImages,
      const std::vector<GPUBuffer>& offImages, int nOffImages,
      const GPUBuffer& CrefImage,
      vector3d *drifts, int nphases, int nx, int ny, int dir,
      float rdistcutoff, float drift_filter_fact)
{
  printf("In determinedrift_2D.\n");
  fflush(stdout);
}

__host__ void separate(int nx, int ny, int nz, int direction, int nphases,
    int norders, std::vector<GPUBuffer>*rawImages, float *sepMatrix)
{
#ifndef NDEBUG
  for (std::vector<GPUBuffer>::iterator i = rawImages->begin();
      i != rawImages->end(); ++i) {
    assert(i->hasNaNs() == false);
  }
#endif
  // Allocate memory for result (have to do this out-of-place)
  std::vector<float*> output(norders * 2 - 1);
  for (std::vector<float*>::iterator i = output.begin(); i != output.end(); ++i) {
    cutilSafeCall(cudaMalloc((void**)&(*i), nz * ny * (nx + 2) *
          sizeof(float)));
  }
  cutilSafeCall(cudaMemcpyToSymbol(const_outputPtrs, &output[0],
        output.size() * sizeof(output[0])));

  // Transfer image pointers in __constant__ array
  std::vector<float*> imgPtrs;
  for (std::vector<GPUBuffer>::iterator i = rawImages->begin(); i != rawImages->end(); ++i) {
    imgPtrs.push_back((float*)i->getPtr());
  }
  cutilSafeCall(cudaMemcpyToSymbol(const_imgPtrs, &imgPtrs[0],
        imgPtrs.size() * sizeof(imgPtrs[0])));

  // Transfer sepMatrix to __constant__ array
  cutilSafeCall(cudaMemcpyToSymbol(const_sepMatrix, &sepMatrix[0],
        (norders * 2 - 1) * nphases * sizeof(sepMatrix[0])));

  // Do the separation step
  int nThreadsX = 16;
  int nThreadsY = 16;
  dim3 nThreads(nThreadsX, nThreadsY, 1);
  int numBlocksX = (int)ceil((float)(nx + 2) / nThreadsX);
  int numBlocksY = (int)ceil((float)ny / nThreadsY);
  dim3 nBlocks(numBlocksX, numBlocksY, 1);
  separate_kernel<<<nBlocks, nThreads>>>( norders, nphases, nx, ny, nz);
  cutilSafeCall(cudaGetLastError());

  // Release the input data pointers and swap the result pointers into
  // the rawImages
  for (int i = 0; i < nphases; ++i) {
    rawImages->at(i).resize(0);
    rawImages->at(i).setPtr((char*)output[i],
        nz * ny * (nx + 2) * sizeof(float), 0);
  }
#ifndef NDEBUG
  for (std::vector<GPUBuffer>::iterator i = rawImages->begin();
      i != rawImages->end(); ++i) {
  //  i->dump(std::cout, nx + 2, 0, nz * ny * (nx + 2)* sizeof(float));
    assert(i->hasNaNs() == false);
  }
#endif
}

__global__ void separate_kernel(int norders, int nphases,
    int nx, int ny, int nz)
{
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  int y = blockIdx.y * blockDim.y  + threadIdx.y;
  int nxy2 = (nx + 2) * ny;
  int offset = y * (nx + 2) + x;
  if (x < nx + 2 && y < ny) {
    for (int i = 0; i < norders * 2 - 1; ++i) {
      float* outBasePtr = const_outputPtrs[i];
      for (int z = 0; z < nz; ++z) {
        float result = 0.0f;
        float* outPtr = outBasePtr + z * nxy2;
        for (int j = 0; j < nphases; ++j) {
          const float* imgBasePtr = const_imgPtrs[j];
          const float* imgPtr = imgBasePtr + z * nxy2;
          float mij= const_sepMatrix[i * nphases + j];
          result +=  mij * imgPtr[offset];
        }
        outPtr[offset] = result;
      }
    }
  }
}

__host__ void makemodeldata(int nx, int ny, int nz, std::vector<GPUBuffer>* bands,
    int norders, vector k0, float dy, float dz,
    std::vector<GPUBuffer>* OTF, short wave, ReconParams *pParams) {
  printf("In makemodeldata.\n");
  fflush(stdout);
}

__host__ void fixdrift_bt_dirs(std::vector<GPUBuffer>* bands, int norders,
    vector3d drift, int nx,int ny, int nz) {
  printf("In fixdrift_bt_dirs.\n");
  fflush(stdout);
}

__host__ void findk0(std::vector<GPUBuffer>* bands, GPUBuffer* overlap0,
    GPUBuffer* overlap1, int nx, int ny, int nz, int norders, vector *k0,
    float dy, float dz, std::vector<GPUBuffer>* OTF, short wave,
    ReconParams * pParams)
{
  int fitorder1;
  int fitorder2;


  fitorder1 = 0;
  if (nz > 1) {
    if (!pParams->bBessel) {
      fitorder2 = 2;
    } else {
      fitorder2 = 1;
    }
  }
  else {
    fitorder2 = 1;
  }

  makeoverlaps(bands, overlap0, overlap1, nx, ny, nz, fitorder1, fitorder2,
      (*k0).x, (*k0).y, dy, dz, OTF, wave, pParams);

  GPUBuffer crosscorr_c(nx * ny * sizeof(cuFloatComplex), 0);
  aTimesConjB(overlap0, overlap1, nx, ny, nz, &crosscorr_c);

  cufftHandle cufftplan;
  int err = cufftPlan2d(&cufftplan, ny, nx, CUFFT_C2C);
  if (CUFFT_SUCCESS != err) {
	printf("cufftPlanxd failed at %s(%d)\n", __FILE__, __LINE__);
    printf("Error code: %d\n", err);
    fflush(stdout);
    exit(-1);
  }
  err = cufftExecC2C(cufftplan, (cuFloatComplex*)crosscorr_c.getPtr(),
      (cuFloatComplex*)crosscorr_c.getPtr(), CUFFT_FORWARD);
  if (CUFFT_SUCCESS != err) {
    printf("cufftExecC2C failed at %s(%d)\n", __FILE__, __LINE__);
    printf("Error code: %d\n", err);
    fflush(stdout);
    exit(-1);
  }
  cufftDestroy(cufftplan);

  GPUBuffer crosscorr(nx * ny * sizeof(float), 0);
  computeIntensities(&crosscorr_c, nx, ny, &crosscorr);
  CPUBuffer intensitiesHost(crosscorr.getSize());
  crosscorr.set(&intensitiesHost, 0, crosscorr.getSize(), 0);
  //  std::cout << "Cross correlation:" << std::endl;
  //  intensitiesHost.dump(std::cout, nx, 0, nx * ny * sizeof(float));

  vector old_k0 = *k0;
  findpeak((float*)intensitiesHost.getPtr(), nx, ny, k0);

  if (old_k0.x < (*k0).x - nx / 2) (*k0).x -= nx;
  if (old_k0.x > (*k0).x + nx / 2) (*k0).x += nx;
  if (old_k0.y < (*k0).y - ny / 2) (*k0).y -= ny;
  if (old_k0.y > (*k0).y + ny / 2) (*k0).y += ny;

  k0->x /= fitorder2;
  k0->y /= fitorder2; /* return k0 of the first order, no matter which fitorder2 is used */
}

__host__ void aTimesConjB(GPUBuffer* overlap0, GPUBuffer* overlap1,
    int nx, int ny, int nz, GPUBuffer* crosscorr_c)
{
  dim3 blockSize;
  blockSize.x = 16;
  blockSize.y = 16;
  blockSize.z = 1;
  dim3 numBlocks;
  numBlocks.x = (int)(ceil((float)nx / blockSize.x));
  numBlocks.y = (int)(ceil((float)ny / blockSize.y));
  numBlocks.z = 1;
  aTimesConjBKernel<<<numBlocks, blockSize>>>(
      (cuFloatComplex*)overlap0->getPtr(), (cuFloatComplex*)overlap1->getPtr(),
      nx, ny, nz,
      (cuFloatComplex*)crosscorr_c->getPtr());
}

__global__ void aTimesConjBKernel(cuFloatComplex* overlap0,
    cuFloatComplex* overlap1, int nx, int ny, int nz,
    cuFloatComplex* crosscorr_c)
{
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int l = blockIdx.y * blockDim.y + threadIdx.y;
  if (k < nx && l < ny) {
    int nxy = nx * ny;
    overlap0 += l * nx + k;
    overlap1 += l * nx + k;
    cuFloatComplex result;
    result.x = 0.0f;
    result.y = 0.0f;
    for (int z = 0; z < nz; ++z) {
      cuFloatComplex Xval = *overlap0;
      cuFloatComplex Yval = *overlap1;
      result.x += Xval.x * Yval.x + Xval.y * Yval.y;
      result.y += -Xval.x * Yval.y + Xval.y * Yval.x;
      overlap0 += nxy;
      overlap1 += nxy;
    }
    crosscorr_c[l * nx + k] = result;
  }
}

__host__ void computeIntensities(GPUBuffer* amplitudes, int nx, int ny,
    GPUBuffer* intensities)
{
  dim3 blockSize;
  blockSize.x = 16;
  blockSize.y = 16;
  blockSize.z = 1;
  dim3 numBlocks;
  numBlocks.x = (int)(ceil((float)nx / blockSize.x));
  numBlocks.y = (int)(ceil((float)ny / blockSize.y));
  numBlocks.z = 1;
  computeIntensitiesKernel<<<numBlocks, blockSize>>>(
      (cuFloatComplex*)amplitudes->getPtr(), nx, ny,
      (float*)intensities->getPtr());
}

__global__ void computeIntensitiesKernel(cuFloatComplex* amplitudes,
    int nx, int ny, float* intensities)
{
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int l = blockIdx.y * blockDim.y + threadIdx.y;
  if (k < nx && l < ny) {
    cuFloatComplex amp = amplitudes[l * nx + k];
    intensities[l * nx + k] = amp.x * amp.x + amp.y * amp.y;
  }
}

__host__ void findpeak(float array[], int sizex, int sizey, vector *peak)
{
  int   xcent=0, ycent=0, i, j;
  float a1, a2, a3, big;

  big = -1e11;
  for(i=0;i<sizey;i++)
    for(j=0;j<sizex;j++)
      if(array[i*sizex+j] > big) {
        big=array[i*sizex+j];
        ycent = i;  xcent = j;
      }

  if(xcent==0)
    a1 = array[ ycent*sizex +  xcent-1+sizex];
  else
    a1 = array[ ycent*sizex +  xcent-1];
  a2 = array[ ycent*sizex +  xcent  ];
  if(xcent==sizex-1)
    a3 = array[ ycent*sizex +  xcent+1-sizex];
  else
    a3 = array[ ycent*sizex +  xcent+1];
  (*peak).x = fitparabola(a1,a2,a3) + xcent;

  if(ycent==0)
    a1 = array[ (ycent-1+sizey)*sizex + xcent ];
  else
    a1 = array[ (ycent-1)*sizex + xcent ];
  a2 = array[ (ycent  )*sizex + xcent ];
  if(ycent==sizey-1)
    a3 = array[ (ycent+1-sizey)*sizex + xcent ];
  else
    a3 = array[ (ycent+1)*sizex + xcent ];
  (*peak).y = fitparabola(a1,a2,a3) + ycent;
}

__host__ float fitparabola( float a1, float a2, float a3 )
{
  float slope,curve,peak;

  slope = 0.5* (a3-a1);         /* the slope at (x=0).  */
  curve = (a3+a1) - 2*a2;       /* (a3-a2)-(a2-a1). The change in slope per unit of x. */
  if( curve == 0 ) {
    printf("no peak: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f\n",a1,a2,a3,slope,curve);
    return( 0.0 );
  }
  peak = -slope/curve;          /* the x value where slope = 0  */
  if( peak>1.5 || peak<-1.5 ) {
    printf("bad peak position: a1=%f, a2=%f, a3=%f, slope=%f, curvature=%f, peak=%f\n",a1,a2,a3,slope,curve,peak);
    return( 0.0 );
  }
  return( peak );
}

__host__ void makeoverlaps(std::vector<GPUBuffer>* bands,
    GPUBuffer* overlap0, GPUBuffer* overlap1, int nx, int ny, int nz,
    int order1, int order2, float k0x, float k0y, float dy, float dz,
    std::vector<GPUBuffer>* OTF, short wave, ReconParams* params)
{
  float order0_2_factor = 1.0f;
  if (nz > 1) {
    order0_2_factor = 5.0f;
    if (params->bBessel)
      order0_2_factor = 4.0f;
  }
  float dkr = 1.0f / (ny * dy);
  float dkz;
  if (dz > 0.0f) {
    dkz = 1.0f / (nz * dz);
  } else {
    dkz = params->dkzotf;
  }
  float krscale = dkr / params->dkrotf;
  float kzscale = dkz / params->dkzotf;
  int rdistcutoff = (int)((params->na * 2.0 / (wave / 1.0e3)) / dkr);
  if (rdistcutoff > nx / 2) {
    rdistcutoff = nx / 2;
  }
  float k0pix = sqrt(k0x * k0x + k0y * k0y);
  float k0mag = k0pix * dkr;
  float lambdaem = (wave / params->nimm) / 1.0e3;
  float lambdaexc = 0.88 * lambdaem;
  float alpha = asin(params->na / params->nimm);
  float beta = asin(k0mag / (2.0 / lambdaexc));
  float betamin = asin(k0mag / (2.0 / lambdaexc) - sin(alpha) *
      SPOTRATIO);
  float zdistcutoff;
  if (!params->bTwolens && !params->bBessel) {
    zdistcutoff = (int)ceil(((1.0 - cos(alpha)) / lambdaem) / dkz);
  }
  else if (params->bBessel) {
    float halfangle;
    float kzExMax;
    kzExMax = 2.0 * params->BesselNA / params->BesselLambdaEx;
    halfangle = acos(k0mag * order2 /(params->norders-1)/ kzExMax);
    zdistcutoff = ceil((kzExMax * sin(halfangle) + (1.0 - cos(alpha)) / lambdaem) / dkz);
  }
  else {
    std::cerr << "Sorry, this program doesn't handel 2-objective mode data\n";
    exit(-1);
  }


  if (zdistcutoff > nz / 2) {
    zdistcutoff = ((nz / 2 - 1) > 0) ? (nz / 2 - 1) : 0;
  }
  printf("order2=%d, rdistcutoff=%d, zdistcutoff=%f\n", order2, rdistcutoff, zdistcutoff);

  float kx = k0x * (order2 - order1);
  float ky = k0y * (order2 - order1);
  float otfcutoff = 0.008;
  if (params->bBessel)
    otfcutoff = 0.01;

  cutilSafeCall(cudaMemset((void*)overlap0->getPtr(), 0,
        nx * ny * nz * sizeof(cuFloatComplex)));
  cutilSafeCall(cudaMemset((void*)overlap1->getPtr(), 0,
        nx * ny * nz * sizeof(cuFloatComplex)));

  // Copy configuration parameters to GPU
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bSuppress_singularities,
        &params->bSuppress_singularities, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_suppression_radius,
        &params->suppression_radius, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bDampenOrder0,
        &params->bDampenOrder0, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bNoKz0,
        &params->bNoKz0, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bFilteroverlaps,
        &params->bFilteroverlaps, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_apodizeoutput,
        &params->apodizeoutput, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bBessel,
        &params->bBessel, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bRadAvgOTF,
        &params->bRadAvgOTF, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_nzotf,
        &params->nzotf, sizeof(int)));
  std::vector<cuFloatComplex*> otfPtrs;

  for (int i = 0; i < params->norders; ++i) {
    otfPtrs.push_back((cuFloatComplex*)(OTF->at(i).getPtr()));
    //    OTF->at(i).dump(std::cout, 100, 0, 100 * sizeof(float));
  }
  cutilSafeCall(cudaMemcpyToSymbol(const_otfPtrs,
        &otfPtrs[0], params->norders * sizeof(cuFloatComplex *),
        0, cudaMemcpyHostToDevice));


  // Set the band ptrs
  cuFloatComplex* band1re;
  cuFloatComplex* band1im;
  cuFloatComplex* band2re;
  cuFloatComplex* band2im;
  //  std::cout << "bands in makeoverlaps:\n";
  if (order1 == 0) {
    //    bands->at(0).dump(std::cout, nx + 2, 0, 2 * (nx + 2) *
    //        sizeof(float));
    band1re = (cuFloatComplex*)bands->at(0).getPtr();
    band1im = 0;
  } else {
    //    bands->at(order1 * 2 - 1).dump(std::cout, nx + 2, 0, 2 * (nx + 2) *
    //        sizeof(float));
    band1re = (cuFloatComplex*)bands->at(order1 * 2 - 1).getPtr();
    //    bands->at(order1 * 2).dump(std::cout, nx + 2, 0, 2 * (nx + 2) *
    //        sizeof(float));
    band1im = (cuFloatComplex*)bands->at(order1 * 2).getPtr();
  }
  // It is assumed that order2 is never 0
  //  bands->at(order2 * 2 - 1).dump(std::cout, nx + 2, 0, 2 * (nx + 2) * sizeof(float));
  band2re = (cuFloatComplex*)bands->at(order2 * 2 - 1).getPtr();
  // bands->at(order2 * 2).dump(std::cout, nx + 2, 0, 2 * (nx + 2) * sizeof(float));
  band2im = (cuFloatComplex*)bands->at(order2 * 2).getPtr();

  // Generate the overlap arrays
  int numThreads = 128;
  dim3 threads(numThreads, 1, 1);
  int numBlocksX = nx / numThreads;
  if (nx % numThreads != 0) {
    ++numBlocksX;
  }
  int numBlocksY = ny;
  int numBlocksZ = 2 * (int)zdistcutoff + 1;
  dim3 blocks(numBlocksX, numBlocksY, numBlocksZ);
  makeOverlaps0Kernel<<<blocks,threads>>>(
      nx, ny, nz, order1, order2, kx, ky, rdistcutoff,
      otfcutoff, zdistcutoff, order0_2_factor, krscale, kzscale,
      band1im, band1re, (cuFloatComplex*)overlap0->getPtr());
  cutilSafeCall(cudaGetLastError());
  makeOverlaps1Kernel<<<blocks,threads>>>(
      nx, ny, nz, order1, order2, kx, ky, rdistcutoff,
      otfcutoff, zdistcutoff, order0_2_factor, krscale, kzscale,
      band2im, band2re, (cuFloatComplex*)overlap1->getPtr());
  cutilSafeCall(cudaGetLastError());

#ifndef NDEBUG
  assert(overlap0->hasNaNs() == false);
  assert(overlap1->hasNaNs() == false);
#endif

  //  std::cout << "Before ffts\n";
  //  std::cout << "overlap 0:\n";
  //  overlap0->dump(std::cout, 2 * nx, 0, 2 * nx * sizeof(float));
  //  std::cout << "overlap 1:\n";
  //  overlap1->dump(std::cout, 2 * nx, 0, 2 * nx * sizeof(float));
  // Do ffts
  cufftResult err;
  cufftHandle cufftplan;
  if (nz > 1) {
    err = cufftPlan3d(&cufftplan, nz, ny, nx, CUFFT_C2C);
  } else {
    err = cufftPlan2d(&cufftplan, ny, nx, CUFFT_C2C);
  }
  if (CUFFT_SUCCESS != err) {
    printf("cufftPlanxd failed\n");
    printf("Error code: %d\n", err);
    fflush(stdout);
    exit(-1);
  }
  err = cufftExecC2C(cufftplan, (cuFloatComplex*)overlap0->getPtr(),
      (cuFloatComplex*)overlap0->getPtr(), CUFFT_INVERSE);
  if (CUFFT_SUCCESS != err) {
    printf("cufftExecC2C failed at %s(%d)\n", __FILE__, __LINE__);
    printf("Error code: %d\n", err);
    fflush(stdout);
    exit(-1);
  }
  err = cufftExecC2C(cufftplan, (cuFloatComplex*)overlap1->getPtr(),
      (cuFloatComplex*)overlap1->getPtr(), CUFFT_INVERSE);
  if (CUFFT_SUCCESS != err) {
    printf("cufftExecC2C failed at %s(%d)\n", __FILE__, __LINE__);
    printf("Error code: %d\n", err);
    fflush(stdout);
    exit(-1);
  }

  cufftDestroy(cufftplan);

#ifndef NDEBUG
  assert(overlap0->hasNaNs() == false);
  assert(overlap1->hasNaNs() == false);
#endif
}

__global__ void makeOverlaps0Kernel(int nx, int ny, int nz,
    int order1, int order2, float kx, float ky,
    float rdistcutoff, float otfcutoff, float zdistcutoff,
    float order0_2_factor, float krscale, float kzscale,
    cuFloatComplex *band1im, cuFloatComplex *band1re,
    cuFloatComplex *overlap0)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j<nx) {
    int i = blockIdx.y;
    int z0 = blockIdx.z - zdistcutoff;
    int x1 = j;
    int y1 = i;
    if (x1 > nx / 2) {
      x1 -= nx;
    }
    if (y1 > ny / 2) {
      y1 -= ny;
    }

    float rdist1 = sqrt((float)x1 * (float)x1 + (float)y1 * (float)y1);
    if (rdist1 <= rdistcutoff) {

      float x12 = x1 - kx;
      float y12 = y1 - ky;
      float rdist12 = sqrt(x12 * x12 + y12 * y12);
      float x21 = x1 + kx;
      float y21 = y1 + ky;
      float rdist21 = sqrt(x21 * x21 + y21 * y21);
      if (rdist12 <= rdistcutoff || rdist21 <= rdistcutoff) { 

        int iin;
        int jin;
        int conj;
        if (j <= nx / 2) {
          iin = i;
          jin = j;
          conj = 0;
        } else {
          jin = nx - j;
          iin = (ny - i) % ny;
          conj = 1;
        }

        if (rdist12 <= rdistcutoff) {
          if (!(z0 == 0 && const_pParams_bNoKz0)) {
            cuFloatComplex otf1 = dev_otfinterpolateMkOvrLps(
                                      const_otfPtrs[order1], x1, y1, krscale, z0, kzscale, 1);
            if (sqrt(otf1.x * otf1.x + otf1.y * otf1.y) > otfcutoff) {
              cuFloatComplex otf12 = dev_otfinterpolateMkOvrLps(
                                         const_otfPtrs[order2], x12, y12, krscale, z0, kzscale, 1);
              if (sqrt(otf12.x * otf12.x + otf12.y * otf12.y) * order0_2_factor > otfcutoff) {
                int z;
                if (conj) {
                  z = -z0;
                } else {
                  z = z0;
                }
                z = (z + nz) % nz;
                int indin = z * (nx / 2 + 1) * ny + iin * (nx / 2 + 1) + jin;
                cuFloatComplex val1re = band1re[indin];
                cuFloatComplex val1im;
                val1im.x = 0.0f;
                val1im.y = 0.0f;
                if (order1 > 0) {
                  val1im = band1im[indin];
                }
                float root = sqrt(otf1.x * otf1.x + otf1.y * otf1.y +
                                  otf12.x * otf12.x + otf12.y * otf12.y);
                cuFloatComplex fact = otf12;
                fact.x /= root;
                fact.y /= root;
                if (conj) {
                  val1re.y *= -1.0;
                  if (order1 > 0) {
                    val1im.y *= -1.0;
                  }
                }
                float temp = val1re.x * fact.x - val1re.y * fact.y;
                val1re.y = val1re.x * fact.y + val1re.y * fact.x;
                val1re.x = temp;
                if (order1 > 0) {
                  temp = val1im.x * fact.x - val1im.y * fact.y;
                  val1im.y = val1im.x * fact.y + val1im.y * fact.x;
                  val1im.x = temp;
                }

                z = (z0 + nz) % nz;
                int indout = z * nx * ny + i * nx + j;
//              if (order1 > 0) {
                overlap0[indout].x = val1re.x - val1im.y;
                overlap0[indout].y = val1re.y + val1im.x;
                // } else {
                //   overlap0[indout] = val1re;
                // }
              }
            }
          }
        }
      }
    }
  }
}

__global__ void makeOverlaps1Kernel(int nx, int ny, int nz,
    int order1, int order2, float kx, float ky,
    float rdistcutoff, float otfcutoff, float zdistcutoff,
    float order0_2_factor, float krscale, float kzscale,
    cuFloatComplex *band2im, cuFloatComplex *band2re,
    cuFloatComplex *overlap1)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j<nx) {
    int i = blockIdx.y;
    int z0 = blockIdx.z - zdistcutoff;
    int x1 = j;
    int y1 = i;
    if (x1 > nx / 2) {
      x1 -= nx;
    }
    if (y1 > ny / 2) {
      y1 -= ny;
    }

    float rdist1 = sqrt((float)x1 * (float)x1 + (float)y1 * (float)y1);
    if (rdist1 < rdistcutoff) {

      float x12 = x1 - kx;
      float y12 = y1 - ky;
      float rdist12 = sqrt(x12 * x12 + y12 * y12);
      float x21 = x1 + kx;
      float y21 = y1 + ky;
      float rdist21 = sqrt(x21 * x21 + y21 * y21);
      if (rdist12 <= rdistcutoff || rdist21 <= rdistcutoff) {

        int iin;
        int jin;
        int conj;
        if (j <= nx / 2) {
          iin = i;
          jin = j;
          conj = 0;
        } else {
          jin = nx - j;
          iin = (ny - i) % ny;
          conj = 1;
        }

        if (rdist21 <= rdistcutoff) {
          if (!(z0 == 0 && const_pParams_bNoKz0)) {
            cuFloatComplex otf2 = dev_otfinterpolateMkOvrLps(
                                       const_otfPtrs[order2], x1, y1, krscale, z0, kzscale, 1);
            if (sqrt(otf2.x * otf2.x + otf2.y * otf2.y) * order0_2_factor > otfcutoff) {
              cuFloatComplex otf21 = dev_otfinterpolateMkOvrLps(
                                          const_otfPtrs[order1], x21, y21, krscale, z0, kzscale, 1);
              if (sqrt(otf21.x * otf21.x + otf21.y * otf21.y) > otfcutoff) {
                int z;
                if (conj) {
                  z = -z0;
                } else {
                  z = z0;
                }
                z = (z + nz) % nz;
                int indin = z * (nx / 2 + 1) * ny + iin * (nx / 2 + 1) + jin;
                cuFloatComplex val2re = band2re[indin];
                cuFloatComplex val2im = band2im[indin];
                float root = sqrt(otf2.x * otf2.x + otf2.y * otf2.y +
                                  otf21.x * otf21.x + otf21.y * otf21.y);
                cuFloatComplex fact = otf21;
                fact.x /= root;
                fact.y /= root;
                if (conj) {
                  val2re.y *= -1.0f;
                  val2im.y *= -1.0f;
                }
                float temp = val2re.x * fact.x - val2re.y * fact.y;
                val2re.y = val2re.x * fact.y + val2re.y * fact.x;
                val2re.x = temp;
                temp = val2im.x * fact.x - val2im.y * fact.y;
                val2im.y = val2im.x * fact.y + val2im.y * fact.x;
                val2im.x = temp;

                z = (z0 + nz) % nz;
                int indout = z * nx * ny + i * nx + j;
                overlap1[indout].x = val2re.x - val2im.y;
                overlap1[indout].y = val2re.y + val2im.x;
              }
            }
          }
        }
      }
    }
  }
}

__device__ cuFloatComplex dev_otfinterpolateMkOvrLps(
    cuFloatComplex * otf, float kx, float ky, float krscale, int kz,
    float kzscale, int mask) {
  cuFloatComplex otfval = make_cuFloatComplex(0.0f, 0.0f);
  if (const_pParams_bRadAvgOTF) {
    int irindex, izindex, indices[2][2];
    float krindex, kzindex;
    float ar, az;

    krindex = sqrt(kx*kx+ky*ky) * krscale;
    kzindex = kz * kzscale;
    if (kzindex<0) kzindex += const_pParams_nzotf;

    irindex = floor(krindex);
    izindex = floor(kzindex);

    ar = krindex - irindex;
    az = kzindex - izindex;
    if (izindex == const_pParams_nzotf-1) {
      indices[0][0] = irindex*const_pParams_nzotf+izindex;
      indices[0][1] = irindex*const_pParams_nzotf;
      indices[1][0] = (irindex+1)*const_pParams_nzotf+izindex;
      indices[1][1] = (irindex+1)*const_pParams_nzotf;
    }
    else {
      indices[0][0] = irindex*const_pParams_nzotf+izindex;
      indices[0][1] = irindex*const_pParams_nzotf+(izindex+1);
      indices[1][0] = (irindex+1)*const_pParams_nzotf+izindex;
      indices[1][1] = (irindex+1)*const_pParams_nzotf+(izindex+1);
    }
    if (mask) {
      otfval.x = (1-ar)*(otf[indices[0][0]].x*(1-az) + otf[indices[0][1]].x*az) +
        ar*(otf[indices[1][0]].x*(1-az) + otf[indices[1][1]].x*az);
      otfval.y = (1-ar)*(otf[indices[0][0]].y*(1-az) + otf[indices[0][1]].y*az) +
        ar*(otf[indices[1][0]].y*(1-az) + otf[indices[1][1]].y*az);
    }
  }
  return otfval;
}

__host__ void fitk0andmodamps(std::vector<GPUBuffer>* bands,
    GPUBuffer* overlap0, GPUBuffer* overlap1, int nx, int ny, int nz,
    int norders, vector *k0, float dy, float dz, std::vector<GPUBuffer>* otf,
    short wave, cuFloatComplex amps[], ReconParams * pParams)
{ 
  int fitorder1 = 0;
  int fitorder2 = 0;
  if (nz > 1) {
	fitorder2 = 2;
    if (pParams->bBessel) 
      fitorder2 = 1;
  }
  else {
    fitorder2 = 1;
  }

  float k0mag = sqrt(k0->x * k0->x + k0->y * k0->y);
  float k0angle = atan2(k0->y, k0->x);

  /* recalculate the overlap arrays at least this first time */
  int redoarrays = (pParams->recalcarrays >= 1);
  float x2 = k0angle;
  cuFloatComplex modamp;
  float amp2 = getmodamp(k0angle, k0mag, bands, overlap0,  overlap1, nx, ny, nz,
      fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);

  /* recalculate the overlap arrays every time only if recalcarrays >= 3*/
  redoarrays = (pParams->recalcarrays >= 3);
  float deltaangle = 0.001;
  float deltamag = 0.1;
  float angle = k0angle + deltaangle;
  float x3 = angle;
  float amp3 = getmodamp(angle, k0mag, bands, overlap0,  overlap1, nx, ny, nz,
      fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);

  float amp1;
  float x1 = 0.0;
  float a;
  if (amp3 > amp2) {
    while(amp3 > amp2) {
      amp1 = amp2;
      x1 = x2;
      amp2 = amp3;
      x2 = x3;
      angle += deltaangle;
      x3 = angle;
      amp3 = getmodamp(angle, k0mag, bands, overlap0, overlap1, nx, ny, nz,
          fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  } else {
    angle = k0angle;
    a = amp3;
    amp3 = amp2;
    amp2 = a;
    a = x3;
    x3 = x2;
    x2 = a;
    while (amp3 > amp2) {
      amp1 = amp2;
      x1 = x2;
      amp2 = amp3;
      x2 = x3;
      angle -= deltaangle;
      x3 = angle;
      amp3 = getmodamp(angle, k0mag, bands, overlap0, overlap1, nx, ny, nz,
          fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  }  /* the maximum of modamp(x) is now between x1 and x3 */
  angle = fitxyparabola(x1, amp1, x2, amp2, x3, amp3);   /* this should be a good angle.  */

  /***** now search for optimum magnitude, at this angle  *****/

  x2 = k0mag;
  amp2 = getmodamp(angle, k0mag, bands, overlap0, overlap1, nx, ny, nz,
      fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);

  float mag = k0mag + deltamag;
  x3 = mag;
  amp3 = getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz,
      fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
  if (amp3 > amp2) {
    while (amp3 > amp2) {
      amp1 = amp2;
      x1 = x2;
      amp2 = amp3;
      x2 = x3;
      mag += deltamag;
      x3 = mag;
      amp3 = getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz,
          fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  } else {
    mag = k0mag;
    a = amp3;
    amp3 = amp2;
    amp2 = a;
    a = x3;
    x3 = x2;
    x2 = a;
    while (amp3 > amp2) {
      amp1 = amp2;
      x1 = x2;
      amp2 = amp3;
      x2 = x3;
      mag -= deltamag;
      x3 = mag;
      amp3 = getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz,
          fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 0);
    }
  }  /* the maximum of modamp(x) is now between x1 and x3 */

  mag = fitxyparabola(x1, amp1, x2, amp2, x3, amp3);  /* this should be a good magnitude.  */

  /* if we were perfectionist we'd iterate for angle again now */

  printf("Optimum modulation amplitude:\n");
  redoarrays = (pParams->recalcarrays>=2);    /* recalculate the d_overlap arrays for optimum modamp fit */
  amp3 = getmodamp(angle, mag, bands, overlap0,  overlap1, nx, ny, nz,
      fitorder1, fitorder2, dy, dz, otf, wave, &modamp, redoarrays, pParams, 1);
  /* one last time, to find the modamp at the optimum k0*/

  float dk = (1/(ny*dy));   /* inverse microns per pixel in data */
  printf("Optimum k0 angle=%f, length=%f, spacing=%f microns\n", angle, mag, 1.0 / (mag * dk));

  k0->x = mag * cos(angle);
  k0->y = mag * sin(angle);
  amps[fitorder2] = modamp;

  /* finally find the modamp for the other orders */
  redoarrays=1;
  if (nz == 1) {
    for (int order = 2; order < norders; ++order) {
      /* assuming that "angle" and "mag" remain the same for every adjacent pair of bands within one direction */
      getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz,
          order - 1, order, dy, dz, otf, wave, &modamp, redoarrays, pParams, 1);
      amps[order] = modamp;
    }
  } else {
    /* 3D */
    for (int order = 1; order < norders; ++order) {
      if (order != fitorder2) {
        getmodamp(angle, mag, bands, overlap0, overlap1, nx, ny, nz,
            0, order, dy, dz, otf, wave, &modamp, redoarrays, pParams, 1);
        amps[order] = modamp;
      }
    }
  }
}

__host__ float fitxyparabola( float x1, float y1, float x2, float y2, float x3, float y3 )
{
  float slope1,slope2,curve,peak,xbar1,xbar2;

  if( x1==x2 || x2==x3 || x3==x1 ) {
    printf("Fit fails; two points are equal: x1=%f, x2=%f, x3=%f\n",x1,x2,x3);
    return( 0.0 );
  }
  xbar1 = 0.5 * (x1 + x2);               /* middle of x1 and x2 */
  xbar2 = 0.5 * (x2 + x3);               /* middle of x2 and x3 */
  slope1 = (y2-y1)/(x2-x1);    /* the slope at (x=xbar1).  */
  slope2 = (y3-y2)/(x3-x2);    /* the slope at (x=xbar2).  */
  curve = (slope2-slope1) / (xbar2-xbar1);       /* The change in slope per unit of x. */
  if( curve == 0 ) {
    printf("Fit fails; no curvature: r1=(%f,%f), r2=(%f,%f), r3=(%f,%f) slope1=%f, slope2=%f, curvature=%f\n",
        x1,y1,x2,y2,x3,y3, slope1,slope2,curve);
    return( 0.0 );
  }

  peak = xbar2 - slope2/curve;          /* the x value where slope = 0  */

  return( peak );
}

__host__ float getmodamp(float kangle, float klength,
    std::vector<GPUBuffer>* bands, GPUBuffer* overlap0, GPUBuffer* overlap1,
    int nx, int ny,int nz, int order1, int order2, float dy, float dz,
    std::vector<GPUBuffer>* otf, short wave, cuFloatComplex* modamp,
    int redoarrays, ReconParams *pParams, int bShowDetail)
{
  vector k1;
  float amp2;
  float corr_coef;
  cuFloatComplex amp_inv;
  cuFloatComplex amp_combo;

  k1.x = klength * cos(kangle);
  k1.y = klength * sin(kangle);
  corr_coef = findrealspacemodamp(bands, overlap0, overlap1, nx, ny, nz, order1, order2, k1, dy, dz, otf,
      wave, modamp, &amp_inv, &amp_combo, redoarrays, pParams);
  amp2 = modamp->x * modamp->x + modamp->y * modamp->y;

  printf(" In getmodamp: angle=%f, mag=%f, amp=%f, phase=%f\n", kangle, klength, sqrt(amp2), get_phase(*modamp));
  if (bShowDetail) {
    printf(" Reverse modamp is: amp=%f, phase=%f\n", 1.0 / cmag(amp_inv), -get_phase(amp_inv));
    printf(" Combined modamp is: amp=%f, phase=%f\n", cmag(amp_combo), get_phase(amp_combo));
    printf(" Correlation coefficient is: %f\n", corr_coef);
  }

  return(amp2);
}

__host__ float findrealspacemodamp(
    std::vector<GPUBuffer>* bands,
    GPUBuffer* overlap0, GPUBuffer* overlap1,
    int nx, int ny, int nz, int order1, int order2,
    vector k0, float dy, float dz,
    std::vector<GPUBuffer>* OTF,
    short wave,
    cuFloatComplex *modamp1, cuFloatComplex *modamp2,
    cuFloatComplex *modamp3, int redoarrays,
    ReconParams *pParams)
{
  if (redoarrays) {
    /* make arrays that contain only the overlapping parts of fourier
       space. Otf-equalize there, set to zero elsewhere  */
    makeoverlaps(bands, overlap0, overlap1, nx, ny, nz, order1, order2,
        k0.x, k0.y, dy, dz, OTF, wave, pParams);
  }

  // Launch reduction kernel
  float kx = k0.x * (order2 - order1);
  float ky = k0.y * (order2 - order1);

  int numRedBlocksX = (int)ceil((float)nx / (float)RED_BLOCK_SIZE_X);
  int numRedBlocksY = (int)ceil((float)ny / (float)RED_BLOCK_SIZE_Y);
  int numRed = numRedBlocksX * numRedBlocksY;
  GPUBuffer XStarY_dev(numRed * sizeof(cuFloatComplex), 0);
  GPUBuffer sumXMag_dev(numRed * sizeof(float), 0);
  GPUBuffer sumYMag_dev(numRed * sizeof(float), 0);
  CPUBuffer XStarY(numRed * sizeof(cuFloatComplex));
  CPUBuffer sumXMag(numRed * sizeof(float));
  CPUBuffer sumYMag(numRed * sizeof(float));

  int numBlocksX = (int)ceil((float)nx / (float)RED_BLOCK_SIZE_X);
  int numBlocksY = (int)ceil((float)ny / (float)RED_BLOCK_SIZE_Y);
  dim3 blocks(numBlocksX, numBlocksY, 1);
  dim3 threads(RED_BLOCK_SIZE_X, RED_BLOCK_SIZE_Y, 1);
  reductionKernel<<<blocks,threads>>>(nx, ny, nz, kx, ky,
      (cuFloatComplex*)overlap0->getPtr(),
      (cuFloatComplex*)overlap1->getPtr(),
      (cuFloatComplex*)XStarY_dev.getPtr(),
      (float*)sumXMag_dev.getPtr(),
      (float*)sumYMag_dev.getPtr());

  // Get partially reduced data
  XStarY_dev.set(&XStarY, 0, XStarY_dev.getSize(), 0);
  sumXMag_dev.set(&sumXMag, 0, sumXMag_dev.getSize(), 0);
  sumYMag_dev.set(&sumYMag, 0, sumYMag_dev.getSize(), 0);

  // Do rest of reduction
  cuFloatComplex XStarYFR = cpuReduce((cuFloatComplex*)XStarY.getPtr(), numRed);
  float sumXMagFR = cpuReduce((float*)sumXMag.getPtr(), numRed);
  float sumYMagFR = cpuReduce((float*)sumYMag.getPtr(), numRed);

  // Compute results
  modamp1->x = XStarYFR.x / sumXMagFR;
  modamp1->y = XStarYFR.y / sumXMagFR;
  modamp2->x = XStarYFR.x / sumYMagFR;
  modamp2->y = -XStarYFR.y / sumYMagFR;
  float tan2beta = 2.0f * sqrt(XStarYFR.x * XStarYFR.x +
      XStarYFR.y * XStarYFR.y) / (sumXMagFR - sumYMagFR);
  float beta = 0.5f * atan(tan2beta);
  if (beta < 0.0f) {
    beta += 0.5f * M_PI;
  }
  float modamp_amp = tan(beta);
  float modamp_arg = atan2(XStarYFR.y, XStarYFR.x);
  modamp3->x = modamp_amp * cos(modamp_arg);
  modamp3->y = modamp_amp * sin(modamp_arg);

  float corr_coef = (XStarYFR.x * XStarYFR.x + XStarYFR.y * XStarYFR.y)
    / (sumXMagFR * sumYMagFR);
  corr_coef = sqrt(corr_coef);

  return corr_coef;
}

__global__ void reductionKernel(
    int nx, int ny, int nz,
    float kx, float ky,
    const cuFloatComplex *overlap0, const cuFloatComplex *overlap1,
    cuFloatComplex *XStarY, float *sumXMag, float *sumYMag) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int strideZ = nx * ny;

  // phase factor
  cuFloatComplex expiphi;
  float angle = 2.0f * M_PI * (
      ((float)i - 0.5f * (float)nx) * kx / (float)nx +
      ((float)j - 0.5f * (float)ny) * ky / (float)ny);
  sincosf(angle, &(expiphi.y), &(expiphi.x));

  // reduction at thread level: each thread sums over z dimension
  cuFloatComplex xsy;
  xsy.x = 0.0f;
  xsy.y = 0.0f;
  float sxm = 0.0f;
  float sym = 0.0f;
  overlap0 += j * nx + i;
  overlap1 += j * nx + i;
  cuFloatComplex czero;
  czero.x = 0.0f;
  czero.y = 0.0f;
  for (int z = 0; z < nz; ++z) {
    cuFloatComplex Xval =
      (i < nx && j < ny) ? overlap0[z * strideZ] : czero;
    cuFloatComplex Yval =
      (i < nx && j < ny) ? overlap1[z * strideZ] : czero;
    xsy.x += Xval.x * Yval.x + Xval.y * Yval.y;
    xsy.y += Xval.x * Yval.y - Xval.y * Yval.x;
    sxm += Xval.x * Xval.x + Xval.y * Xval.y;
    sym += Yval.x * Yval.x + Yval.y * Yval.y;
  }
  float temp = xsy.x * expiphi.x - xsy.y * expiphi.y;
  xsy.y = xsy.x * expiphi.y + xsy.y * expiphi.x;
  xsy.x = temp;

  // reduction at thread block level: summation over tile in x-y plane
  // this is a pretty naive implementation.
  __shared__ cuFloatComplex xsyShrd[RED_BLOCK_SIZE_Y * RED_BLOCK_SIZE_X];
  __shared__ float sxmShrd[RED_BLOCK_SIZE_Y * RED_BLOCK_SIZE_X]; 
  __shared__ float symShrd[RED_BLOCK_SIZE_Y * RED_BLOCK_SIZE_X]; 
  int ltid = threadIdx.y * RED_BLOCK_SIZE_X + threadIdx.x; 
  xsyShrd[ltid] = xsy;
  sxmShrd[ltid] = sxm; 
  symShrd[ltid] = sym; 
  __syncthreads();
  int s;
  for (s = RED_BLOCK_SIZE_Y * RED_BLOCK_SIZE_X / 2; s > 0; s >>= 1) {
    if (ltid < s) {
      xsyShrd[ltid].x += xsyShrd[ltid + s].x;
      xsyShrd[ltid].y += xsyShrd[ltid + s].y;
      sxmShrd[ltid] += sxmShrd[ltid + s];
      symShrd[ltid] += symShrd[ltid + s];
    }
    __syncthreads();
  }

  // write result back
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    int blockIndex = blockIdx.y * gridDim.x + blockIdx.x;
    XStarY[blockIndex] = xsyShrd[0];
    sumXMag[blockIndex] = sxmShrd[0];
    sumYMag[blockIndex] = symShrd[0];
  }
}

__host__ void filterbands(int dir, std::vector<GPUBuffer>* bands,
    const std::vector<vector>& k0, int ndirs, int norders,
    std::vector<GPUBuffer>& otf, float dy, float dz,
    const std::vector<std::vector<cuFloatComplex> >& amp,
    const std::vector<float>& noiseVarFactors, int nx, int ny, int nz,
    short wave, ReconParams* pParams)
{
  float *ampmag2;
  cuFloatComplex *conjamp;
  float rdistcutoff;
  int order, order2, dir2, *zdistcutoff;
  float apocutoff, zapocutoff;
  float dkr, dkz, krscale, kzscale, k0mag, k0pix;
  /* int iin, jin, conj, xyind, ind, z0, iz, z; */
  float lambdaem, lambdaexc, alpha, beta, betamin, wiener;

  wiener = pParams->wiener*pParams->wiener;
  dkr = (1/(ny*dy));   /* inverse microns per pixel in data */
  if (dz>0)
    dkz = (1/(nz*dz));   /* inverse microns per pixel in data */
  else
    dkz = pParams->dkzotf;
  krscale = dkr / pParams->dkrotf;   /* ratio of radial direction pixel scales of data and otf */
  kzscale = dkz / pParams->dkzotf;   /* ratio of axial direction pixel scales of data and otf */
  k0pix =  sqrt(k0[0].x*k0[0].x + k0[0].y*k0[0].y);   /* k0 magnitude (for highest order) in pixels */
  k0mag = k0pix * dkr;   /* k0 magnitude (for highest order) in inverse microns */
  lambdaem = (wave/pParams->nimm)/1000.0;  /* emission wavelength in the sample, in microns */
  lambdaexc = 0.88* lambdaem;;  /* 0.88 approximates a typical lambdaexc/lambdaem  */
  alpha = asin(pParams->na/pParams->nimm);  /* aperture angle of objectives */
  beta = asin(k0mag/(2/lambdaexc));   /* angle of center of side illumination beams */
  betamin = asin((k0mag/(2/lambdaexc)) -sin(alpha)*SPOTRATIO);   /* angle of inner edge of side illumination beams */
  rdistcutoff = (pParams->na*2/(wave/1000.0)) / dkr;    /* OTF support radial limit in data pixels */
  if (rdistcutoff>nx/2) rdistcutoff=nx/2;

  /* 080201: zdistcutoff[0] depends on options -- single or double lenses */
  zdistcutoff = (int *) malloc(norders * sizeof(int));
  if (!pParams->bTwolens && !pParams->bBessel) {
    zdistcutoff[0] = (int) ceil(((1-cos(alpha))/lambdaem) / dkz);    /* OTF support axial limit in data pixels */
    zdistcutoff[norders-1] = 1.3*zdistcutoff[0];    /* approx max axial support limit of the OTF of the high frequency side band */
    if (norders>=3)
      for (order=1;order<norders-1;order++)
        zdistcutoff[order] = (1+lambdaem/lambdaexc)*zdistcutoff[0];       /* axial support limit of the OTF of the medium frequency side band(s?) */
  }
  else if (pParams->bBessel) {
    float kzExMax, halfangle;
    kzExMax = 2 *pParams->BesselNA / pParams->BesselLambdaEx;

    zdistcutoff[0] = (int) rint((kzExMax + (1-cos(alpha))/lambdaem) / dkz);    /* OTF support axial limit in data pixels */
    printf("norders=%d, zdistcutoff[%d]=%d\n", norders, 0 ,zdistcutoff[0]);
    for (order=1; order<norders; order++) {
      halfangle = acos(k0mag * order / (norders-1) / kzExMax);
      zdistcutoff[order] = ceil((kzExMax * sin(halfangle) + (1.0 - cos(alpha)) / lambdaem) / dkz);
      printf("zdistcutoff[%d]=%d\n", order ,zdistcutoff[order]);
    }
  }
  else {  /* two lenses */
    zdistcutoff[0] = (int) ceil(1.02*(2/lambdaem + 2/lambdaexc) / dkz);  /* 1.02 is just a safety margin */
    zdistcutoff[norders-1] = (int) ceil(1.02*(2/lambdaem + 2*cos(beta)/lambdaexc) / dkz);    /* approx max axial support limit of the OTF of the high frequency side band */
    if (norders==3) {
      zdistcutoff[1] =  (int) ceil(1.02*(2/lambdaem + (1+cos(betamin))/lambdaexc) / dkz); /* axial support limit of the OTF of the medium frequency side band */
    }
    else if (norders>3)
      for (order=1;order<norders-1;order++) {
        float a;
        a = ((float)order)/(norders-1);
        zdistcutoff[order] = 1.1*((1-a)*zdistcutoff[0] + a*zdistcutoff[norders-1]);       /* axial support limit of the OTF of the medium frequency side bands */ /* 1.1 is a blur margin */
      }
  }

  for (order=0;order<norders;order++) {
    if (zdistcutoff[order]>=nz/2) zdistcutoff[order]=((nz/2-1) > 0 ? (nz/2-1) : 0);
    /* printf("order=%d, rdistcutoff=%f, zdistcutoff=%d\n", order, rdistcutoff, zdistcutoff[order]); */
  }

  apocutoff = rdistcutoff+ k0pix*(norders-1);

  if (pParams->bTwolens || pParams->bBessel)
    zapocutoff = zdistcutoff[0];
  else
    zapocutoff = zdistcutoff[1];

  ampmag2 = (float *) malloc(norders * sizeof(float));
  conjamp = (cuFloatComplex *) malloc(norders * sizeof(cuFloatComplex));
  for (order=0;order<norders;order++) {
    ampmag2[order] = amp[dir][order].x * amp[dir][order].x +
      amp[dir][order].y * amp[dir][order].y;
    conjamp[order] = amp[dir][order];
    conjamp[order].y *= -1;
  }

  /////////////////////////////////////
  // Move Data into Constant memory
  /////////////////////////////////////
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bSuppress_singularities,
        &pParams->bSuppress_singularities, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_suppression_radius,
        &pParams->suppression_radius, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bDampenOrder0,
        &pParams->bDampenOrder0, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bNoKz0, &pParams->bNoKz0,
        sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bFilteroverlaps,
        &pParams->bFilteroverlaps, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_apodizeoutput,
        &pParams->apodizeoutput, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_apoGamma,
        &pParams->apoGamma, sizeof(float)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bBessel,
        &pParams->bBessel, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_bRadAvgOTF,
        &pParams->bRadAvgOTF, sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_pParams_nzotf, &pParams->nzotf,
        sizeof(int)));
  cutilSafeCall(cudaMemcpyToSymbol(const_wiener, &wiener,
        sizeof(float)));

  cutilSafeCall(cudaMemcpyToSymbol(const_zdistcutoff, zdistcutoff,
        norders*sizeof(int), 0, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaMemcpyToSymbol(const_ampmag2, ampmag2,
        norders * sizeof(float), 0, cudaMemcpyHostToDevice));

  // Explicitly calculate mag2 of amp for all orders
  float * ampmag2_alldirs;
  ampmag2_alldirs = (float *) malloc(ndirs*norders*sizeof(float));
  for (dir2=0; dir2<ndirs; dir2++) {
    for (order2=0;order2<norders;order2++) {
      ampmag2_alldirs[dir2*norders+order2] =
        amp[dir2][order2].x * amp[dir2][order2].x +amp[dir2][order2].y * amp[dir2][order2].y;
    }
  }
  cutilSafeCall(cudaMemcpyToSymbol(const_ampmag2_alldirs,
        ampmag2_alldirs, ndirs* norders * sizeof(float), 0,
        cudaMemcpyHostToDevice));
  free(ampmag2_alldirs);

  cutilSafeCall(cudaMemcpyToSymbol(const_conjamp, conjamp,
        norders * sizeof(cuFloatComplex), 0, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaMemcpyToSymbol(const_noiseVarFactors,
        &noiseVarFactors[0], ndirs* norders * sizeof(float),
        0, cudaMemcpyHostToDevice));
  cutilSafeCall(cudaMemcpyToSymbol(const_k0, &k0[0],
        ndirs * sizeof(float2),
        0, cudaMemcpyHostToDevice));
  //  DM 13/12/2012: The OTFs look fine at this point.
  //  dumpBands(&otf, 128, 257, 1);
  //  exit(0);
  std::vector<cuFloatComplex*> otfPtrs;
  for (std::vector<GPUBuffer>::iterator i = otf.begin();
      i != otf.end(); ++i) {
    otfPtrs.push_back((cuFloatComplex*)i->getPtr());
  }
  cutilSafeCall(cudaMemcpyToSymbol(const_otfPtrs, &otfPtrs[0], 
        norders*sizeof(cuFloatComplex *), 0, cudaMemcpyHostToDevice));

#ifndef NDEBUG
  ///////////////////////////////////////////////////////
  // Allocate space for temporary buffers
  ///////////////////////////////////////////////////////
  cuFloatComplex *dev_tempbandplus; // cuFloatComplex array of size nx*ny*nz
  cuFloatComplex *dev_scale;        // cuFloatComplex array of size nx*ny*nz
  cutilSafeCall(cudaMalloc((void**) &dev_tempbandplus,
        nx*ny*nz*sizeof(cuFloatComplex)));
  cutilSafeCall(cudaMemset((void*) dev_tempbandplus, 0,
        nx*ny*nz*sizeof(cuFloatComplex)));

  cutilSafeCall(cudaMalloc((void**) &dev_scale, nx*ny*nz*sizeof(cuFloatComplex)));
  cutilSafeCall(cudaMemset((void*) dev_scale, 0, nx*ny*nz*sizeof(cuFloatComplex)));
#endif

  cuFloatComplex * dev_bandptr, * dev_bandptr2;
  for (order=0;order<norders;order++) {
    if (order==0) {
      dev_bandptr = (cuFloatComplex*)bands->at(0).getPtr();
    }
    else {
      dev_bandptr = (cuFloatComplex*)bands->at(2*order-1).getPtr();     /* bands contains only the data of one direction -- dir*/
      dev_bandptr2 = (cuFloatComplex*)bands->at(2*order).getPtr();
    }


    //
    // KERNEL 1
    //
    int nThreads = 128;
    int NZblock = 2*zdistcutoff[order]+1;
    int NYblock = ny;
    int NXblock = (int) ceil( (float)nx/2./nThreads );
    dim3 grid(NXblock, NYblock, NZblock);
    dim3 block(nThreads, 1, 1);

    filterbands_kernel1<<<grid,block>>>(dir, ndirs, order, norders, nx, ny, nz,
        rdistcutoff, zapocutoff, apocutoff, krscale, kzscale,
        dev_bandptr, dev_bandptr2, false);
    cutilSafeCall(cudaGetLastError());

    filterbands_kernel1<<<grid,block>>>(dir, ndirs, order, norders, nx, ny, nz,
        rdistcutoff, zapocutoff, apocutoff, krscale, kzscale,
        dev_bandptr, dev_bandptr2, true);
    cutilSafeCall(cudaGetLastError());

#ifndef NDEBUG
    CPUBuffer dev_scale_debug(nx * ny * nz * sizeof(cuFloatComplex));
    cudaMemcpy(dev_scale_debug.getPtr(), dev_scale, nx * ny
        * nz * sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    assert(dev_scale_debug.hasNaNs() == false);
#endif

    //
    // KERNEL 3
    //
    if ((nz-zdistcutoff[order]) > (zdistcutoff[order]+1)) {
      NZblock = (nz-zdistcutoff[order]) - (zdistcutoff[order]+1);
      NXblock = (int) ceil( (float)(nx+2)/2./nThreads );
      dim3 grid2(NXblock, NYblock, NZblock);
      filterbands_kernel3<<<grid2,block>>>(order, nx, ny, nz,
          dev_bandptr, dev_bandptr2);
      cutilSafeCall(cudaGetLastError());
    }

  } /* for order */

  free(zdistcutoff);
  free(ampmag2);
  free(conjamp);
  return;
}

__global__ void filterbands_kernel1(int dir, int ndirs, int order, int norders, int nx, int ny, 
    int nz, float rdistcutoff, float zapocutoff, float apocutoff, float krscale, float kzscale,
    cuFloatComplex * dev_bandptr, cuFloatComplex * dev_bandptr2, bool bSecondEntry)
{

  float kx, ky, rdist1, rdistabs, apofact;
  kx = order * const_k0[dir].x;
  ky = order * const_k0[dir].y;

  // compute x1, y1, z0 based on block and thread indices
  int x1 = blockIdx.x * blockDim.x + threadIdx.x + 1;
  if (x1 <= nx/2) {

    if (bSecondEntry)
      x1 -= nx/2;
    int y1 = blockIdx.y - ny/2;
    int z0 = blockIdx.z - const_zdistcutoff[order];

    float xabs, yabs;
    int iin, jin, conj, xyind, ind, iz, z;
    cuFloatComplex scale, bandreval, bandimval, bandplusval, bandminusval;

    cuFloatComplex otf1, otf2;
    float weight, sumweight, dampfact;
  
    /*x1, y1 are coords within each band to be scaled */
    if (x1 >= 0) {
      /* integer coords of actual arrays to be filtered */
      iin = y1; 
      jin = x1;
      conj = 0;
    } else {
      iin = -y1;
      jin = -x1;
      conj = 1;
    }

    if (order==0 && conj) return;   /* For center band only the actual pixels need to be filtered */

    if (iin<0) {
      iin += ny;
    }
    xyind = iin * (nx / 2 + 1) + jin;

    rdist1 = sqrtf(x1*x1+y1*y1);  /* dist from center of band to be filtered */

    // The following variable masks the whole computation
    if (rdist1<=rdistcutoff) { /* is x1,y1 within the theoretical lateral OTF support of 
                                  the data that is to be scaled? */
      xabs=x1+kx;   /* (floating point) coords rel. to absolute fourier space, with */
      yabs=y1+ky;   /* the absolute origin=(0,0) after the band is shifted by k0 */
      rdistabs = sqrt(xabs*xabs + yabs*yabs);  // used later for apodization calculation
      otf1 = dev_otfinterpolate(const_otfPtrs[order], x1, y1, krscale, z0, kzscale);

      weight = otf1.x * otf1.x + otf1.y * otf1.y;
      if (order!= 0) weight *= const_ampmag2[order];
      dampfact = 1. / const_noiseVarFactors[dir*norders+order];
    
      // this one is thread dependent ... from the rdist calculation
      if (const_pParams_bSuppress_singularities && order != 0 && rdist1 <=const_pParams_suppression_radius)
        dampfact *= dev_suppress(rdist1);
    
      // these next two are not thread dependent
      else if (!const_pParams_bDampenOrder0 && const_pParams_bSuppress_singularities && order ==0)
        dampfact *= dev_suppress(rdist1);
    
      else if (const_pParams_bDampenOrder0 && order ==0)
        dampfact *= dev_order0damping(rdist1, z0, rdistcutoff, const_zdistcutoff[0]);
    
      // if no kz=0 plane is used:
      if (order==0 && z0==0 && const_pParams_bNoKz0) dampfact = 0;
    
      weight *= dampfact;
      sumweight=weight;

      int dir2, order2;
      float amp2mag2, kx2, ky2, rdist2;
      float x2, y2;
      for (dir2=0; dir2<ndirs; dir2++) {
        for (order2=-(norders-1); order2<norders; order2++) {
          if (dir2==dir && order2==order) continue;
          if (!const_pParams_bFilteroverlaps && !(order2==0 && order==0)) continue; /* bFilteroverlaps is always true except when (during debug) generating an unfiltered exploded view */
          amp2mag2 = const_ampmag2_alldirs[dir2*norders+abs(order2)];
          kx2 = order2 * const_k0[dir2].x;
          ky2 = order2 * const_k0[dir2].y;
          x2 = xabs-kx2; /* coords rel to shifted center of band 2 */
          y2 = yabs-ky2;
          rdist2 = sqrt(x2*x2+y2*y2);       /* dist from center of band 2 */
    
          if (rdist2<rdistcutoff) {
      
            otf2 = dev_otfinterpolate(const_otfPtrs[abs(order2)], x2, y2, krscale, z0, kzscale);
            weight = dev_mag2(otf2) / const_noiseVarFactors[dir2*norders+abs(order2)];
            if (order2 != 0) weight *= amp2mag2;
      
            if (const_pParams_bSuppress_singularities && order2 != 0 && rdist2 <= const_pParams_suppression_radius)
              weight *= dev_suppress(rdist2);
      
            else if (!const_pParams_bDampenOrder0 && const_pParams_bSuppress_singularities && order2 ==0)
              weight *= dev_suppress(rdist2);
      
            else if (const_pParams_bDampenOrder0 && order2==0)
              weight *= dev_order0damping(rdist2, z0, rdistcutoff, const_zdistcutoff[0]);
      
            if (const_pParams_bNoKz0 && order2==0 && z0==0) weight = 0.0f;
      
            sumweight += weight;
          }
        }
      }
  
      sumweight += const_wiener;
      scale.x = dampfact *   otf1.x/sumweight;
      scale.y = dampfact * (-otf1.y)/sumweight;

      if (const_pParams_apodizeoutput) {
        float rho, zdistabs;
        zdistabs = abs(z0);

        if (zapocutoff > 0) {  /* 3D case */
          if (!const_pParams_bBessel)
			rho = sqrt((rdistabs / apocutoff) * (rdistabs / apocutoff) +
					   (zdistabs / zapocutoff) * (zdistabs / zapocutoff));
		  else {
            float rhox, rhoy, rhoz;
            rhox = xabs/rdistcutoff; //apocutoff * 1.6f;
            rhoy = yabs/apocutoff;
            rhoz = zdistabs/zapocutoff;
            rho = sqrt(rhox*rhox + rhoy*rhoy + rhoz*rhoz);
          }
        }
        else         /* 2D case */
          rho = sqrt((rdistabs/apocutoff)*(rdistabs/apocutoff));

        if (rho > 1.f) rho = 1.0f;

        if (const_pParams_apodizeoutput == 1)    /* cosine-apodize */
          apofact = cos((M_PI*0.5f)* rho);
        else if (const_pParams_apodizeoutput == 2)
          apofact = 1.0f - rho;
        // apofact = __powf(1.0f - rho, const_pParams_apoGamma);
        scale.x *= apofact;
        scale.y *= apofact;
      }
      /* What we want is to use mag2 for the weights, as you have done, and
       * then set  scale = conjugate(otf1)/sumweight */

      /* separate (for this pixel) the even and odd "bands" into the true
       * plus and minus bands */
      /* apply the scale to the plus band only (this will apply it to the
       * minus band as well by symmetry?) */
      /* reassemble into the even and odd "bands" */
      if (conj) {
        z = -z0;
      } else {
        z = z0;
      }
      /* coords of the fourier space arrays have the origin of fourier space
       * at (0,0,0) */
      iz = (z + nz) % nz;
      /* index of the corresponding point in the input array */
      ind = iz*((nx/2+1)*ny) + xyind;
      if (order == 0) {
        // if (!conj) // Condition "order == 0 && conj" has been ruled out earlier
        dev_bandptr[ind] = cuCmulf(dev_bandptr[ind], scale);
      }
      else {
        scale = cuCmulf(scale, const_conjamp[order]); /* not invamp: the 1/|amp| factor is
                                                         taken care of by including ampmag2 in the weights */
        bandreval = dev_bandptr[ind];
        bandimval = dev_bandptr2[ind];
        if (conj) {
          bandreval.y *= -1.0f;
          bandimval.y *= -1.0f;
        }
        /* bandplus = bandre + i bandim */
        bandplusval.x = bandreval.x - bandimval.y;
        bandplusval.y = bandreval.y + bandimval.x;
        /* bandminus = bandre - i bandim */
        bandminusval.x = bandreval.x + bandimval.y;
        bandminusval.y = bandreval.y - bandimval.x;
        /* scale only the bandplus part - bandminus will take care 
           of itself because of symmetry (?) */
        bandplusval = cuCmulf(bandplusval, scale);

        bandreval.x = 0.5f*( bandplusval.x + bandminusval.x);
        bandreval.y = 0.5f*( bandplusval.y + bandminusval.y);
        bandimval.x = 0.5f*( bandplusval.y - bandminusval.y);
        bandimval.y = 0.5f*(-bandplusval.x + bandminusval.x);
        if (conj) {
          bandreval.y *= -1.f;
          bandimval.y *= -1.f;
        }
        dev_bandptr[ind] = bandreval;
        dev_bandptr2[ind] = bandimval;
      }
    }
    else { //if (rdisk1>...)
      iz = (z0 + nz) % nz;
      ind = iz * ((nx / 2 + 1) * ny) + xyind;
      dev_bandptr[ind] = make_cuFloatComplex(0.f, 0.f);
      if (order != 0)
        dev_bandptr2[ind] = make_cuFloatComplex(0.f, 0.f);
    }
  }
  return;
}



__global__ void filterbands_kernel3(int order, int nx, int ny, int nz,
    cuFloatComplex * dev_bandptr, cuFloatComplex * dev_bandptr2) {
//! Clear everything above and below zdistcutoff to 0

  int x1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (x1 < nx/2+1) {
    int y1 = blockIdx.y;
    int z0 = blockIdx.z + const_zdistcutoff[order] + 1;

    int ind = z0*((nx/2+1)*ny) + y1 * (nx/2+1) + x1;
    dev_bandptr[ind] = make_cuFloatComplex(0.f, 0.f);
    if (order !=0)
      dev_bandptr2[ind] = make_cuFloatComplex(0.f, 0.f);
  }
  return;
}

__device__ cuFloatComplex dev_otfinterpolate(cuFloatComplex * otf, float kx, float ky, float krscale, int kz, float kzscale)
  /* (kx, ky, kz) is Fourier space coords with origin at kx=ky=kz=0 and going  betwen -nx(or ny,nz)/2 and +nx(or ny,nz)/2 */
{
  cuFloatComplex otfval = make_cuFloatComplex(0.f, 0.f);
  // This should be rewritten using Textures for the interpolation. It will be much easier and faster!
  if (const_pParams_bRadAvgOTF) {
    int irindex, izindex, indices[2][2];
    float krindex, kzindex;
    float ar, az;

    krindex = sqrt(kx*kx+ky*ky) * krscale;
    kzindex = kz * kzscale;
    if (kzindex<0) kzindex += const_pParams_nzotf;

    irindex = floor(krindex);
    izindex = floor(kzindex);

    ar = krindex - irindex;
    az = kzindex - izindex;  // az is always 0 for 2D case, and it'll just become a 1D interp

    if (izindex == const_pParams_nzotf-1) {
      indices[0][0] = irindex*const_pParams_nzotf+izindex;
      indices[0][1] = irindex*const_pParams_nzotf;
      indices[1][0] = (irindex+1)*const_pParams_nzotf+izindex;
      indices[1][1] = (irindex+1)*const_pParams_nzotf;
    }
    else {
      indices[0][0] = irindex*const_pParams_nzotf+izindex;
      indices[0][1] = irindex*const_pParams_nzotf+(izindex+1);
      indices[1][0] = (irindex+1)*const_pParams_nzotf+izindex;
      indices[1][1] = (irindex+1)*const_pParams_nzotf+(izindex+1);
    }
    otfval.x = (1-ar)*(otf[indices[0][0]].x*(1-az) + otf[indices[0][1]].x*az) +
      ar*(otf[indices[1][0]].x*(1-az) + otf[indices[1][1]].x*az);
    otfval.y = (1-ar)*(otf[indices[0][0]].y*(1-az) + otf[indices[0][1]].y*az) +
      ar*(otf[indices[1][0]].y*(1-az) + otf[indices[1][1]].y*az);
  }
  return otfval;
}

__device__ float dev_order0damping(float radius, float zindex, int rlimit, int zlimit)
{
  float rfraction, zfraction;

  rfraction = radius/rlimit;
  zfraction = fabs(zindex/zlimit);

  return rfraction*rfraction + zfraction*zfraction*zfraction;
}

__device__ float dev_mag2(cuFloatComplex x)
{
  return x.x*x.x+x.y*x.y;
}

__device__ float dev_suppress(float x)
{
  float x6,out;
  x6 = x*x*x;
  x6 *= x6;
  out = 1.0/(1+20000/(x6+20));
  return out;
}

__host__ void assemblerealspacebands(int dir, GPUBuffer* outbuffer,
    GPUBuffer* bigbuffer, std::vector<GPUBuffer>* bands, int ndirs,
    int norders, const std::vector<vector>& k0, int nx, int ny, int nz,
    float zoomfact, int z_zoom, float expfact)
{
  float fact,* dev_coslookup,* dev_sinlookup;
  int order;

  /* Allocate temporaries */    
  cutilSafeCall(cudaMalloc((void **) &dev_coslookup,
                           (int)(rint(nx*zoomfact)*rint(ny*zoomfact)*sizeof(float))));
  cutilSafeCall(cudaMalloc((void **) &dev_sinlookup,
                           (int)(rint(nx*zoomfact)*rint(ny*zoomfact)*sizeof(float))));

  fact = expfact/0.5;  /* expfact is used for "exploded view".  For normal reconstruction expfact = 1.0  */

  int nThreads = 128;
  int NZblock = nz;
  int NYblock = ny;
  int NXblock = nx/nThreads;
  if (nx%nThreads) NXblock ++;

  dim3 grid(NXblock, NYblock, NZblock);
  dim3 block(nThreads, 1, 1);

  printf("moving centerband\n");
  cutilSafeCall(cudaMemset((void*) bigbuffer->getPtr(), 0,
        unsigned (rint(zoomfact*nx)*rint(zoomfact*ny)*(z_zoom*nz)*sizeof(cuFloatComplex))));
  move_kernel<<<grid,block>>>(
      (cuFloatComplex*)bands->at(0).getPtr(),
      (cuFloatComplex*)bands->at(0).getPtr(),
      0, (cuFloatComplex*)bigbuffer->getPtr(), nx, ny, nz, zoomfact, z_zoom);

  cufftHandle myGPUPlan;
  cufftResult cuFFTErr = cufftPlan3d(&myGPUPlan, (int) (z_zoom*nz), (int) rint(zoomfact*ny),
                                     (int) rint(zoomfact*nx), CUFFT_C2C);
  if (cuFFTErr!=CUFFT_SUCCESS) {
    if (cuFFTErr == CUFFT_ALLOC_FAILED)
      printf("\n*** In assemblerealspacebands(), CUFFT failed to allocate GPU or CPU memory\n");
    throw std::runtime_error("CUFFT plan creation failed");
  }

  /* transform it */
  printf("re-transforming centerband\n");
  cuFFTErr = CUFFT_SUCCESS;
  cuFFTErr = cufftExecC2C(myGPUPlan,
      (cuFloatComplex*)bigbuffer->getPtr(),
      (cuFloatComplex*)bigbuffer->getPtr(),
      CUFFT_INVERSE);
  if (cuFFTErr!=CUFFT_SUCCESS) printf("Error in cufftExecC2C: %d\n", cuFFTErr);

  printf("inserting centerband\n");
  NZblock = (int)(z_zoom*nz);
  NYblock = (int) rint(zoomfact*ny);
  NXblock = (int) ceil((zoomfact*nx)/nThreads);
  dim3 grid2(NXblock, NYblock, NZblock);
  write_outbuffer_kernel1<<<grid2,block>>>((cuFloatComplex*)bigbuffer->getPtr(),
      (float*)outbuffer->getPtr(), (int) rint(zoomfact*nx));

  printf("centerband assembly completed\n");

  for (order=1; order < norders; order ++) {
    float k0x, k0y;
    /* move side bands to bigbuffer, fill in with zeroes */
    printf("moving order %d\n",order); 
    cutilSafeCall(cudaMemset((void*) bigbuffer->getPtr(), 0,   unsigned (rint(zoomfact*nx)*rint(zoomfact*ny)*(z_zoom*nz)*sizeof(cuFloatComplex))));
    move_kernel<<<grid,block>>>((cuFloatComplex*)bands->at(2*order-1).getPtr(),
        (cuFloatComplex*)bands->at(2*order).getPtr(), order,
        (cuFloatComplex*)bigbuffer->getPtr(), nx, ny, nz, zoomfact, z_zoom);

    // transform it into real space
    cuFFTErr = cufftExecC2C(myGPUPlan, (cuFloatComplex*)bigbuffer->getPtr(),
        (cuFloatComplex*)bigbuffer->getPtr(), CUFFT_INVERSE);
    if (cuFFTErr!=CUFFT_SUCCESS) printf("Error in cufftExecC2C: %d\n", cuFFTErr);

    /***** For 3D, prepare 2D array of sines and cosines first, then loop over z. ******/
    k0x = k0[dir].x*((float)order);
    k0y = k0[dir].y*((float)order);

    NZblock = 1;
    NYblock = (int)rint (zoomfact*ny);
    NXblock = (int) ceil(zoomfact*nx/nThreads);
    dim3 grid3(NXblock, NYblock, NZblock);
    cos_sin_kernel<<<grid3,block>>>(k0x,  k0y, fact, dev_coslookup, dev_sinlookup, (int)(zoomfact*nx));
    write_outbuffer_kernel2<<<grid2, block>>>(dev_coslookup,
        dev_sinlookup, (cuFloatComplex*)bigbuffer->getPtr(), 
        (float*)outbuffer->getPtr(), (int) (zoomfact*nx));

    printf("order %d sideband assembly completed\n", order);
  } /* for (order =...) */

  /* Free memory */
  cutilSafeCall(cudaFree((void *) dev_coslookup));
  cutilSafeCall(cudaFree((void *) dev_sinlookup));
  cufftDestroy(myGPUPlan);
  return;
}

__global__ void move_kernel(cuFloatComplex *inarray1, cuFloatComplex *inarray2, int order, 
    cuFloatComplex *outarray, int nx, int ny, int nz, float zoomfact, int z_zoom)
{
  int     xdim, ydim, zdim, nxy, nxyout;

  xdim=rint(zoomfact*nx); ydim=rint(zoomfact*ny);
  zdim = z_zoom*nz;
  nxy = (nx/2+1)*ny;
  nxyout = xdim*ydim;

  // compute x1, y1, z0 based on block and thread indices
  int x = blockIdx.x * blockDim.x + threadIdx.x;
  if (x<nx) {
    x -= nx/2-1;
    int y = blockIdx.y - (ny/2-1);
    int z = (nz>1) ? blockIdx.z - (nz/2-1) : 0;

    int  indin, indout, xout, yout, zout, conj;
    cuFloatComplex valre, valim, val;

    xout = x;     /* xout,yot,zout = (non-centered) output coords with zoomed-up dims and origin of fourier space at (0,0,0) */
    if (xout<0) xout += xdim;
    yout = y;
    if (yout<0) yout += ydim;
    zout = z;
    if (zout<0) zout += zdim;
    indout = zout*nxyout + yout*xdim + xout;

    if (x<0) {    /* now xyz get turned into coords of the (half fourier space) input arrays */
      x = -x;
      y = -y;
      z = -z;
      conj = 1;
    }
    else
      conj = 0;

    if (y<0) y += ny;
    if (z<0) z += nz;
    indin = z*nxy + y*(nx/2+1) + x;

    if (order == 0) {
      val = inarray1[indin];
      if (conj)
        val.y *= -1;
    }
    else {
      valre = inarray1[indin];
      valim = inarray2[indin];
      if (conj) {
        valre.y *= -1;
        valim.y *= -1;
      }
      val.x = valre.x - valim.y;
      val.y = valre.y + valim.x;
    }

    outarray[indout] = val;
  }
  return;
}
__global__ void write_outbuffer_kernel1(cuFloatComplex * bigbuffer, float * outbuffer, int nx) {

  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j<nx) {
    int i = blockIdx.y;
    int k = blockIdx.z;
    int NXlocal = nx; // why was "gridDim.x*blockDim.x" used?
    int NYlocal = gridDim.y;
    int ind = k*NXlocal*NYlocal + i*NXlocal + j;
    outbuffer[ind] += bigbuffer[ind].x;
  }
}

__global__ void write_outbuffer_kernel2(float * coslookup, float * sinlookup, 
      cuFloatComplex * bigbuffer, float * outbuffer, int nx) {

  int j = blockIdx.x * blockDim.x + threadIdx.x;
  if (j<nx) {
    int i = blockIdx.y;
    int k = blockIdx.z;
    int NXlocal = nx; // why was "gridDim.x*blockDim.x" used?
    int NYlocal = gridDim.y;
    int indxy = i*NXlocal + j;
    int ind = k*NXlocal*NYlocal + indxy;
    outbuffer[ind] += (bigbuffer[ind].x * 2.0*coslookup[indxy] - bigbuffer[ind].y * 2.0*sinlookup[indxy]);
  }
}

__global__ void cos_sin_kernel(float k0x, float k0y, float fact, float * coslookup, float * sinlookup, int nx) {

  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y;

  if (j<nx) {
    int ind = i*(gridDim.x*blockDim.x) + j;
    int NXlocal = nx; // why was "gridDim.x*blockDim.x" used?
    int NYlocal = gridDim.y;
    float angle = fact * M_PI * ((j-NXlocal/2)*k0x/NXlocal + (i-NYlocal/2)*k0y/NYlocal);
    coslookup[ind] = cos(angle);
    sinlookup[ind] = sin(angle);
  }
}

__host__ void computeAminAmax(const GPUBuffer* data, int nx, int ny, int nz,
    float* min, float* max)
{
  int numElems = nx * ny * nz;

  int blockSize = 1024;
  int numBlocks = 100;
  GPUBuffer maxPartialResult(numBlocks * sizeof(float), 0);
  GPUBuffer minPartialResult(numBlocks * sizeof(float), 0);
  computeAminAmax_kernel<<<numBlocks, blockSize,
    2 * blockSize * sizeof(float)>>>((const float*)data->getPtr(),
        numElems,
        (float*)maxPartialResult.getPtr(),
        (float*)minPartialResult.getPtr());
  CPUBuffer maxPartialResultHost(maxPartialResult.getSize());
  CPUBuffer minPartialResultHost(minPartialResult.getSize());
  maxPartialResult.set(&maxPartialResultHost, 0,
      maxPartialResult.getSize(), 0);
  minPartialResult.set(&minPartialResultHost, 0,
      minPartialResult.getSize(), 0);
  const float* maxArray = (const float*)maxPartialResultHost.getPtr();
  const float* minArray = (const float*)minPartialResultHost.getPtr();
  *max = -10000.0;
  *min = 10000.0;
  for (int i = 0; i < numBlocks; ++i) {
    if (maxArray[i] > *max) {
      *max = maxArray[i];
    }
    if (minArray[i] < *min) {
      *min = minArray[i];
    }
  }
}

__global__ void computeAminAmax_kernel(const float* data, int numElems,
    float* maxPartialResult, float* minPartialResult)
{
  volatile extern __shared__ float s_maxmin[];
  s_maxmin[threadIdx.x] = -10000.0f;
  s_maxmin[threadIdx.x + blockDim.x] = 10000.0f;
  for (int i = threadIdx.x;
      i < (int)ceil((float)numElems / (blockDim.x * gridDim.x));
      i += (blockDim.x * gridDim.x)) {
    if (i < numElems) {
      float d = data[i];
      if (d > s_maxmin[threadIdx.x]) {
        s_maxmin[threadIdx.x] = d;
      }
      if (d < s_maxmin[threadIdx.x + blockDim.x]) {
        s_maxmin[threadIdx.x + blockDim.x] = d;
      }
    }
  }
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      if (s_maxmin[threadIdx.x + s] > s_maxmin[threadIdx.x]) {
        s_maxmin[threadIdx.x] = s_maxmin[threadIdx.x + s];
      }
      if (s_maxmin[blockDim.x + threadIdx.x + s] < s_maxmin[blockDim.x + threadIdx.x]) {
        s_maxmin[blockDim.x + threadIdx.x] = s_maxmin[blockDim.x + threadIdx.x + s];
      }
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    maxPartialResult[blockIdx.x] = s_maxmin[0];
    minPartialResult[blockIdx.x] = s_maxmin[blockDim.x];
  }
}


// compute the mean above the background using GPU reduction
__host__ double meanAboveBackground_GPU(GPUBuffer &img, int nx, int ny, int nz)
{
  unsigned nThreads = 1024;
  unsigned nBlocks = (unsigned) ceil( nx*ny*nz /(float) nThreads/2 );
  unsigned smemSize = nThreads * sizeof(double);

  // used for holding intermediate reduction results; one for each thread block
  GPUBuffer d_intres(nBlocks * sizeof(double), 0);

  summation_kernel<<<nBlocks, nThreads, smemSize>>>((float *) img.getPtr(),
                                                    (double *) d_intres.getPtr(), nx*ny*nz);
  // download intermediate results to host:
  CPUBuffer intRes(d_intres);
  double sum=0;
  double *p=(double *)intRes.getPtr();
  for (int i=0; i<nBlocks; i++)
    sum += *p++;

  float mean = sum/((nx-2)*ny*nz);

  GPUBuffer d_counter(nBlocks * sizeof(unsigned), 0);
  smemSize = nThreads * (sizeof(double) + sizeof(unsigned));
  sumAboveThresh_kernel<<<nBlocks, nThreads, smemSize>>>((float *) img.getPtr(),
                                                         (double *) d_intres.getPtr(),
                                                         (unsigned *) d_counter.getPtr(),
                                                         mean, nx*ny*nz);
  
  // download intermediate results to host:
  CPUBuffer counter(d_counter);
  intRes = d_intres;
  sum=0;
  unsigned count = 0;
  p=(double *)intRes.getPtr();
  unsigned *pc = (unsigned *) counter.getPtr();
  for (int i=0; i<nBlocks; i++) {
    sum += *p++;
    count += *pc++;
  }

  // printf("mean=%e, sum=%e, count=%d\n", mean, sum, count);
  return sum/count;
}

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
// (Copied from reduction_kernel.cu of CUDA samples)
template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

__global__ void summation_kernel(float * img, double * intRes, int n)
// Copied from CUDA "reduction" sample code reduce4()
{
  double *sdata = SharedMemory<double>();

  unsigned tid = threadIdx.x;
  unsigned ind = blockIdx.x * blockDim.x*2 + threadIdx.x;

  double mySum= (ind < n) ? img[ind] : 0;

  if (ind + blockDim.x < n)
    mySum += img[ind + blockDim.x];

  sdata[tid] = mySum;
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=blockDim.x/2; s>32; s>>=1) {
    if (tid < s) {
      sdata[tid] += sdata[tid + s];
    }
    __syncthreads();
  }

  if (tid < 32) {
    // now that we are using warp-synchronous programming (below)
    // we need to declare our shared memory volatile so that the compiler
    // doesn't reorder stores to it and induce incorrect behavior.
    volatile double *smem = sdata;

    // Assuming blockSize is > 64:
    smem[tid] += smem[(tid + 32)];
    smem[tid] += smem[(tid + 16)];
    smem[tid] += smem[(tid +  8)];
    smem[tid] += smem[(tid +  4)];
    smem[tid] += smem[(tid +  2)];
    smem[tid] += smem[(tid +  1)];
  }
  // write result for this block to global mem
  if (tid == 0) intRes[blockIdx.x] = sdata[0];
}


__global__ void sumAboveThresh_kernel(float * img, double * intRes, unsigned * counter, float thresh, int n)
// Adapted from CUDA "reduction" sample code reduce4()
{
// Size of shared memory allocated is nThreads * (sizeof(double) + sizeof(unsigned))
// The first nThreads * sizeof(double) bytes are used for image intensity sum;
// the next nThreads * sizeof(unsigned) bytes are for counting pixels whose intensity is > thresh
  double *sdata = SharedMemory<double>();
  unsigned *count = (unsigned *) (sdata + blockDim.x);

  unsigned tid = threadIdx.x;
  unsigned ind = blockIdx.x * blockDim.x*2 + threadIdx.x;

  double mySum= 0;
  unsigned myCount = 0;
  if (ind < n && img[ind] > thresh) {
    mySum = img[ind]; 
    myCount ++;
  }

  unsigned ind2 = ind + blockDim.x;
  if (ind2 < n && img[ind2] > thresh) {
    mySum += img[ind2];
    myCount ++;
  }

  sdata[tid] = mySum;
  count[tid] = myCount;
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=blockDim.x/2; s>32; s>>=1) {
    if (tid < s) {
      sdata[tid] += sdata[tid + s];
      count[tid] += count[tid + s];
    }
    __syncthreads();
  }

  if (tid < 32) {
    volatile double *smem = sdata;
    volatile unsigned *cmem = count;

    smem[tid] += smem[(tid + 32)];
    smem[tid] += smem[(tid + 16)];
    smem[tid] += smem[(tid +  8)];
    smem[tid] += smem[(tid +  4)];
    smem[tid] += smem[(tid +  2)];
    smem[tid] += smem[(tid +  1)];
    cmem[tid] += cmem[(tid + 32)];
    cmem[tid] += cmem[(tid + 16)];
    cmem[tid] += cmem[(tid +  8)];
    cmem[tid] += cmem[(tid +  4)];
    cmem[tid] += cmem[(tid +  2)];
    cmem[tid] += cmem[(tid +  1)];
  }
  // write result for this block to global mem
  if (tid == 0) {
    intRes[blockIdx.x] = sdata[0];
    counter[blockIdx.x] = count[0];
  }
}

__host__ void rescale_GPU(GPUBuffer &img, int nx, int ny, int nz, float scale)
{
  unsigned nThreads = 1024;
  unsigned nBlocks = (unsigned) ceil( nx*ny*nz / (float) nThreads );
  scale_kernel<<<nBlocks, nThreads>>>((float *) img.getPtr(), scale, nx*ny*nz);
#ifndef NDEBUG
  std::cout<< "rescale_GPU(): " << cudaGetErrorString(cudaGetLastError()) << std::endl;
#endif
}

__global__ void scale_kernel(float * img, double factor, int n)
{
  unsigned ind = blockIdx.x * blockDim.x + threadIdx.x;
  if (ind < n)
    img[ind] *= factor;
}

// // // // //! Add 2 to the X dimensions of "*in" and "*out" 
// // // // __global__ void deskew_kernel(float *in, int nx, int ny, int nz,
// // // //                               float *out, int nxOut, int extraShift,
// // // //                               double deskewFactor, float fillVal)
// // // // {
// // // //   unsigned xout = blockIdx.x * blockDim.x + threadIdx.x;
// // // //   unsigned yout = blockIdx.y;
// // // //   unsigned zout = blockIdx.z;

// // // //   if (xout < nxOut) {
// // // //     float xin = (xout - nxOut/2.+extraShift) - deskewFactor*(blockIdx.z-nz/2.) + nx/2.;

// // // //     unsigned indout = zout * (nxOut+2) * ny + yout * (nxOut+2) + xout;
// // // //     if (xin >= 0 && xin < nx-1) {

// // // //       // 09-03-2013 Very important lesson learned:
// // // //       // the (unsigned int) casting has be placed right there because
// // // //       // otherwise, the entire express would evaluate as floating point and
// // // //       // there're only 24-bit mantissa, so any odd index that's > 16777216 would
// // // //       // inaccurately rounded up. int or unsigned does not have the 24-bit limit.
// // // //       unsigned indin = zout * (nx+2) * ny + yout * (nx+2) + (unsigned int) floor(xin);

// // // //       float offset = xin - floor(xin);
// // // //       out[indout] = (1-offset) * in[indin] + offset * in[indin+1];
// // // //     }
// // // //     else
// // // //       out[indout] = fillVal;
// // // //   }
// // // // }

// // // // __host__ void deskew_GPU(std::vector<GPUBuffer> * pImgs, int nx, int ny, int nz,
// // // // 						 float deskewAngle, float dz_prior_to, 
// // // //                          float dr, int extraShift, float fillVal)
// // // // {
// // // //   assert(ny >= nx);
// // // //   int newNx = ny;

// // // //   if (deskewAngle <0) deskewAngle += 180.;
// // // //   float deskewFactor = cos(deskewAngle * M_PI/180.) * dz_prior_to / dr;   // cos() or sin() ?? --lin

// // // //   GPUBuffer outBuf((newNx+2) * ny * nz * sizeof(float), 0);
// // // //   outBuf.setToZero();

// // // //   dim3 block(512, 1, 1);
// // // //   unsigned nxBlocks = (unsigned) ceil(newNx / (float) block.x);
// // // //   dim3 grid(nxBlocks, ny, nz);

// // // //   for (std::vector<GPUBuffer>::iterator inBuf=pImgs->begin(); inBuf != pImgs->end(); inBuf++) {
// // // //     assert(inBuf->hasNaNs() == false);
// // // // 	deskew_kernel<<<grid, block>>>((float *) inBuf->getPtr(), nx, ny, nz,
// // // // 								   (float *) outBuf.getPtr(), newNx,
// // // // 								   extraShift, deskewFactor, fillVal);
// // // // 	inBuf->resize((newNx+2) * ny * nz * sizeof(float));
// // // // 	outBuf.set(&(*inBuf), 0, outBuf.getSize(), 0);
// // // // #ifndef NDEBUG
// // // // 	std::cout<< "deskew_GPU(): " << cudaGetErrorString(cudaGetLastError()) << std::endl;
// // // //     assert(outBuf.hasNaNs() == false);
// // // // #endif
// // // //   }
// // // // }

