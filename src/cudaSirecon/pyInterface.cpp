#include "cudaSireconImpl.h"  // must come first!
#include "SIM_reconstructor.hpp"

extern "C" {
SIM_Reconstructor* SR_new(int nx, int ny, int nImages,
                          const char* configFileName) {
  // std::string cfg(configFileName);
  return new SIM_Reconstructor(nx, ny, nImages, configFileName);
}

void SR_setRaw(SIM_Reconstructor* sr, const float* const raw_data, int nx,
               int ny, int nz) {
  CImg<float> raw_image(raw_data, nx, ny, nz);
  sr->setRaw(raw_image, 0, 0);
  sr->loadAndRescaleImage(0, 0);
  sr->setCurTimeIdx(0);
  sr->processOneVolume();
}

void SR_getResult(SIM_Reconstructor* sr, float* result) {
  return sr->getResult(result);
}

ReconParams& SR_getReconParams(SIM_Reconstructor* sr) {
  return sr->getReconParams();
}
ImageParams& SR_getImageParams(SIM_Reconstructor* sr) {
  return sr->getImageParams();
}
}
