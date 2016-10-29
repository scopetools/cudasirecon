# CUDA_SIMrecon
Lin Shao's 3-beam SIM reconstruction software


#######################################################################
# Build instruction for CUDA SIMrecon project
# 1. Prerequisites
#  1.1 Cmake (> 2.8.11)
#  1.2 CUDA SDK (>5.0)
#  1.3a. If Linux or Mac OS X, GCC (< 4.8)
#  1.3b. If Windows, Visual C++ (> 2012)
#  1.4 Boost libraries (> 1.48)
#  1.5 IVE libraries (download from http://msg.ucsf.edu/IVE/Download/index.html). Place the priism-4.2.9/ folder under this folder where this file is located
#
# 2. Make a build dir (assumed to be "build" under the main folder where this 
#    CMakeLists.txt is located) and cd into it
#
# 3. On Linux or Mac at a shell prompt, type:
#  $ cmake -D CMAKE_BUILD_TYPE=Release ..
#  On Windows, one extra flag is needed. Type the following in a 
#  VS build environment:
#  > cmake -D CMAKE_BUILD_TYPE=Release -G "NMake Makefiles" ..
#  Make sure there isn't any error or warning messages. 
#  Always delete everything in the build dir before re-run cmake
#
# 4. Type "make" on Linux/Mac or "nmake" on Window to build the executable
#
# 5. If building is successful, an executable cudaSireconDriveris generated in
#  build/cudaSirecon. Run the test to make sure it works:
#  $ ../testData/job.sh
#
#######################################################################
