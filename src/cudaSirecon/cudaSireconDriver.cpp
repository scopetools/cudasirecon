#include "cudaSireconImpl.h"
#include "SIM_reconstructor.hpp"


int main(int argc, char **argv)
{
  try {
    SIM_Reconstructor myreconstructor(argc, argv);

    for (int it = 0; it < myreconstructor.getNTimes(); ++it) {
      for (int iw = 0; iw < 1; ++iw) {
        myreconstructor.setFile(it, iw);
        myreconstructor.loadAndRescaleImage(it, iw);
        myreconstructor.setCurTimeIdx(it);
        if (myreconstructor.processOneVolume())
          myreconstructor.writeResult(it, iw);
      }
    }

    myreconstructor.closeFiles();
  }
  catch (std::exception &e) {
    std::cerr << "\n!!Error occurred: " << e.what() << std::endl;
    return 0;
  }
  return 0;
}
