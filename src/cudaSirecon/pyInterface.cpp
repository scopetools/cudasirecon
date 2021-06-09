#include "cudaSireconImpl.h"  // must come first!
#include "SIM_reconstructor.hpp"

extern "C" {

SIM_Reconstructor* SR_new_from_shape(int nx, int ny, int nImages,
                          const char* configFileName) {
  return new SIM_Reconstructor(nx, ny, nImages, configFileName);
}

void SR_setRaw(SIM_Reconstructor* sr,
               const float* const raw_data,
               int nx, int ny, int nz) {
  CImg<float> raw_image(raw_data, nx, ny, nz);
  sr->setRaw(raw_image);
}
void SR_loadAndRescaleImage(SIM_Reconstructor* sr, int it, int iw) {
  sr->loadAndRescaleImage(it, iw);
}
void SR_setCurTimeIdx(SIM_Reconstructor* sr, int it) {
  sr->setCurTimeIdx(it);
}
void SR_processOneVolume(SIM_Reconstructor* sr) {
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

// these functions only used to emulate cli

SIM_Reconstructor* SR_new_from_argv(int argc, char **argv) {
  return new SIM_Reconstructor(argc, argv);
}

void SR_setFile(SIM_Reconstructor* sr, int it, int iw) {
  sr->setFile(it, iw);
}
void SR_writeResult(SIM_Reconstructor* sr, int it, int iw) {
  sr->writeResult(it, iw);
}
void SR_closeFiles(SIM_Reconstructor* sr) {
  sr->closeFiles();
}

}
