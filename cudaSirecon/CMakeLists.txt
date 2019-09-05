include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  ${CUDA_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers"
  )

set(CMAKE_MACOSX_RPATH 1)

if(WIN32)# OR APPLE)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  include_directories(  ${Boost_INCLUDE_DIRS} )
endif()

if(WIN32)
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/INCLUDE")
	include_directories(${PRIISM_INCLUDE_PATH})
    link_directories ( ${PRIISM_LIB_PATH} )
    if(WITH_TIFF)
        include_directories( ${TIFF_INCLUDE_DIR} )
        link_directories ( ${TIFF_LIBRARY} )
    endif(WITH_TIFF)
elseif(APPLE)
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/darwin64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/darwin64/INCLUDE")
	include_directories(${PRIISM_INCLUDE_PATH})
else()
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/linux64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/linux64/INCLUDE")
	include_directories(${PRIISM_INCLUDE_PATH})
endif(WIN32)

message(STATUS "PRIISM_LIB_PATH: " ${PRIISM_LIB_PATH})

find_library(IMLIB imlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
find_library(IVELIB ive PATHS ${PRIISM_LIB_PATH} REQUIRED)


if(WIN32)
  find_library(IMLIB libimlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
  find_library(IVELIB libive PATHS ${PRIISM_LIB_PATH} REQUIRED)

  if ( ${MSVC_VERSION} GREATER 1800)
    find_library(LEG_STDIO legacy_stdio_definitions)
    message(STATUS "adding legacy stdio definitions for >VS14.0 " ${LEG_STDIO})
  endif()

  find_library(BLAS openblas libblas 
    PATHS 
    "${CONDA_PREFIX}/Library/lib"
    "${CMAKE_SOURCE_DIR}/lapack/win64"
    REQUIRED
  )
  find_library(LAPACK liblapack
    PATHS
    "${CONDA_PREFIX}/Library/lib"
    "${CMAKE_SOURCE_DIR}/lapack/win64"
    REQUIRED
  )
else()
  find_library(IMLIB imlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
  find_library(IVELIB ive PATHS ${PRIISM_LIB_PATH} REQUIRED)
endif()

message(STATUS "IMLIB: " ${IMLIB})
message(STATUS "IVELIB: " ${IVELIB})
message(STATUS "LEG_STDIO: " ${LEG_STDIO})
message(STATUS "BLAS: " ${BLAS})
message(STATUS "LAPACK: " ${LAPACK})

find_package(CUDA)
message(STATUS "CUDA: " ${CUDA_CUFFT_LIBRARIES})

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

  if ( ${MSVC_VERSION} GREATER 1800)
    target_link_libraries(
      cudaSireconDriver
      cudaSirecon
      gpuFunctions
      Buffer
      ${LEG_STDIO}
      ${IMLIB}
      ${IVELIB}
      ${LAPACK}
      ${BLAS}
    )
  else()
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
  endif()
else(WIN32)
  target_link_libraries(
    cudaSirecon
    gpuFunctions
    Buffer
    ${IMLIB}
    ${IVELIB}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_REGEX_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    blas
    lapack
    ${TIFF_LIBRARIES}
    X11
  )
  if(NOT APPLE)
    target_link_libraries(
      cudaSirecon
      rt
    )
  endif()
  target_link_libraries(
    cudaSireconDriver
    cudaSirecon
  )
endif(WIN32)
CUDA_ADD_CUFFT_TO_TARGET(cudaSirecon)
CUDA_ADD_CUFFT_TO_TARGET(cudaSireconDriver)

# added for make install to work in conda
set(HEADERS
  cudaSirecon.h
  cudaSireconImpl.h
  gpuFunctions.h
)

install(
  TARGETS cudaSireconDriver cudaSirecon
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(
  FILES ${HEADERS} DESTINATION include
)

# install(
#   FILES ${IMLIB} ${IVELIB} DESTINATION lib
# )
