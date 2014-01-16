#include "cudaSireconImpl.h"
#include "SIM_reconstructor.hpp"

int main(int argc, char **argv)
{
  SIM_Reconstructor myreconstructor(argc, argv);

  for (int it = 0; it < myreconstructor.getNTimes(); ++it) {
    for (int iw = 0; iw < 1; ++iw) {
      myreconstructor.loadAndRescaleImage(it, iw);
      myreconstructor.setCurTimeIdx(it);
      myreconstructor.processOneVolume();
      myreconstructor.writeResult(it, iw);
    }
  }

#ifndef __SIRECON_USE_TIFF__
  myreconstructor.closeFiles();
#endif
  return 0;
}
