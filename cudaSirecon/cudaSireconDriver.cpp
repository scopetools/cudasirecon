#include "cudaSireconImpl.h"
#include "SIM_reconstructor.hpp"

#if defined(_MSC_VER) && (_MSC_VER >= 1900)
// for IMLIB with vs >2015
// https://stackoverflow.com/questions/30412951/unresolved-external-symbol-imp-fprintf-and-imp-iob-func-sdl2
FILE _iob[] = {*stdin, *stdout, *stderr};
extern "C" FILE * __cdecl __iob_func(void)
{
    return _iob;
}
#endif

int main(int argc, char **argv)
{
  try {
    SIM_Reconstructor myreconstructor(argc, argv);

    for (int it = 0; it < myreconstructor.getNTimes(); ++it) {
      for (int iw = 0; iw < 1; ++iw) {
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
