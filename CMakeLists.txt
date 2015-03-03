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


cmake_minimum_required (VERSION 2.8)

project (cudaSIMRecon)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

ENABLE_TESTING()

set(VERSION_MAJOR "0")
set(VERSION_MINOR "0")
set(VERSION_PATCH "1")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

######################################################################
#
# Set permissions before adding subdirectories
#
######################################################################
set(SCI_GROUP_WRITE GROUP_WRITE)
set(SCI_WORLD_FILE_PERMS WORLD_READ)
set(SCI_WORLD_PROGRAM_PERMS WORLD_READ WORLD_EXECUTE)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (APPLE)
# Apple's default clang complier does not support OpenMP
#  set(CMAKE_C_COMPILER gcc)
#  set(CMAKE_CXX_COMPILER g++)
#  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -undefined dynamic_lookup -flat_namespace")
#  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup -flat_namespace")
endif()


# Building the documentation
option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
if(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR 
      "Doxygen is needed to build the documentation. Please
      install it correctly")
    endif()
    #-- Configure the Template Doxyfile for our specific project
    configure_file(Doxyfile.in
      ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)
    #-- Add a custom target to run Doxygen when ever the project is built
    add_custom_target (Docs ALL 
      COMMAND ${DOXYGEN_EXECUTABLE}
      ${PROJECT_BINARY_DIR}/Doxyfile
      SOURCES
      ${PROJECT_BINARY_DIR}/Doxyfile)
endif()

find_package(CUDA)

if (WITH_TIFF)
  if (WIN32)
    set(TIFF_INCLUDE_DIR "c:/libtiff")
    set(TIFF_LIBRARY "c:/libtiff")
  endif(WIN32)
  find_package(TIFF)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -D__SIRECON_USE_TIFF__")
endif()

if(WIN32)
  find_package( Boost )
else()
  find_package( Boost REQUIRED COMPONENTS program_options filesystem regex system)
endif(WIN32)

if(WIN32 OR APPLE)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  include_directories(  ${Boost_INCLUDE_DIRS} )
endif()

if(APPLE)
  find_package(X11)
  link_directories ( ${X11_INCLUDE_DIR}/../lib )
endif()

set(CUDA_NVCC_FLAGS
  -gencode=arch=compute_20,code=sm_20;-gencode=arch=compute_30,code=sm_30;-gencode=arch=compute_35,code=sm_35)

if(WIN32)
  set(CUDA_NVCC_FLAGS
     ${CUDA_NVCC_FLAGS};-I${Boost_INCLUDE_DIRS};--use-local-env;--cl-version=2013)
endif(WIN32)

add_subdirectory(Buffers)
add_subdirectory(cudaSirecon)