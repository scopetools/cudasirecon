include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  ${CUDA_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers"
  )


if(WIN32)# OR APPLE)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  include_directories(  ${Boost_INCLUDE_DIRS} )
endif()

if(WIN32)
	set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/LIB")
	include_directories("${CMAKE_SOURCE_DIR}/IVE/win64/INCLUDE/")
    link_directories ( ${PRIISM_LIB_PATH} )
    if(WITH_TIFF)
        include_directories( ${TIFF_INCLUDE_DIR} )
        link_directories ( ${TIFF_LIBRARY} )
    endif(WITH_TIFF)
elseif(APPLE)
## Assuming the builder's shell environment is already set up to use 64-bit Priism
	set(PRIISM_LIB_PATH "$ENV{IVE_BASE}/Darwin64/LIB")
	include_directories("$ENV{IVE_BASE}/Darwin64/INCLUDE/")
else()
## Assuming the builder's shell environment is already set up to use 64-bit Priism
	set(PRIISM_LIB_PATH "$ENV{IVE_BASE}/Linux/x86_64/LIB")
	include_directories("$ENV{IVE_BASE}/Linux/x86_64/INCLUDE/")
endif(WIN32)


if(WIN32)
  message(STATUS "PRIISM_LIB_PATH: " ${PRIISM_LIB_PATH})
  # find_library(IMLIB libimlib.a
  #   PATHS ${PRIISM_LIB_PATH}
  #   REQUIRED)
  # find_library(IVELIB ive
  #   PATHS ${PRIISM_LIB_PATH}
  #   REQUIRED)
else(WIN32)
  find_library(IMLIB imlib
    PATHS ${PRIISM_LIB_PATH}
    REQUIRED)
  find_library(IVELIB ive
    PATHS ${PRIISM_LIB_PATH}
    REQUIRED)
 endif(WIN32)

# if (APPLE)
###############################################
#### Lin: The path to the cuda libraries is hardcoded here. Let me know
#### if this works for you.
###############################################
#   set(CUDA_LIBRARIES "/Developer/NVIDIA/CUDA-5.0/lib/libtlshook.dylib" ${CUDA_LIBRARIES} )
# endif ()



CUDA_ADD_LIBRARY(
  gpuFunctions
#  SHARED
  gpuFunctionsImpl.cu
  )

if (NOT APPLE)
set (SHARED_OR_NOT SHARED)
endif()

if (WITH_TIFF AND WIN32)
CUDA_ADD_LIBRARY(
  cudaSirecon
  ${SHARED_OR_NOT}
  cudaSirecon.cpp
  boostfs.cpp
  tiffhandle.cpp
  )
elseif(APPLE)
  CUDA_ADD_LIBRARY(
    cudaSirecon
    SHARED
    cudaSirecon.cpp
  )
else ()
  CUDA_ADD_LIBRARY(
    cudaSirecon
#    ${SHARED_OR_NOT}
    cudaSirecon.cpp
  )
endif()

CUDA_ADD_EXECUTABLE(
  cudaSireconDriver
  cudaSireconDriver.cpp
)

add_dependencies(cudaSireconDriver
  cudaSirecon
  gpuFunctions
  Buffer
  )
if(WIN32)
  find_library(BLAS LIBBLAS PATHS "${CMAKE_SOURCE_DIR}/lapack/win64" REQUIRED)
  find_library(LAPACK LIBLAPACK PATHS "${CMAKE_SOURCE_DIR}/lapack/win64" REQUIRED)
  if(WITH_TIFF)
    target_link_libraries(
      cudaSirecon
      gpuFunctions
      Buffer
      ${LAPACK}
      ${BLAS}
      libtiff
    )
  else(WITH_TIFF)
    target_link_libraries(
      cudaSirecon
      gpuFunctions
      Buffer
      # ${IMLIB}
      # ${IVELIB}
	  libimlib
	  libive
      ${LAPACK}
      ${BLAS}
    )
  endif(WITH_TIFF)

  target_link_libraries(
    cudaSireconDriver
    cudaSirecon
    gpuFunctions
    Buffer
    ${IMLIB}
    ${IVELIB}
    ${LAPACK}
    ${BLAS}
  )
else(WIN32)
  target_link_libraries(
    cudaSirecon
    gpuFunctions
    Buffer
    ${IMLIB}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    blas
    lapack
    ${TIFF_LIBRARIES}
    X11
#    nvidia-ml
  )
  target_link_libraries(
    cudaSireconDriver
    cudaSirecon
  )
endif(WIN32)
CUDA_ADD_CUFFT_TO_TARGET(cudaSirecon)
CUDA_ADD_CUFFT_TO_TARGET(cudaSireconDriver)
