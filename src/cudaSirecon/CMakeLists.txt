include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}"
  ${CUDA_INCLUDE_DIRS}
  "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers"
)

if(WIN32)
  link_directories ( ${Boost_LIBRARY_DIRS} )
  link_directories ( ${PRIISM_LIB_PATH} )
  link_directories ( ${TIFF_LIBRARY} )
  include_directories(  ${Boost_INCLUDE_DIRS} )
  include_directories( ${TIFF_INCLUDE_DIR} )
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/INCLUDE")
else()
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/linux64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/linux64/INCLUDE")
endif(WIN32)

include_directories(${PRIISM_INCLUDE_PATH})


message(STATUS "PRIISM_LIB_PATH: " ${PRIISM_LIB_PATH})
find_library(IMLIB NAMES imlib libimlib PATHS ${PRIISM_LIB_PATH})
find_library(IVELIB NAMES ive libive PATHS ${PRIISM_LIB_PATH})

message(STATUS "IMLIB: " ${IMLIB})
message(STATUS "IVELIB: " ${IVELIB})
message(STATUS "LEG_STDIO: " ${LEG_STDIO})

find_library(LAPACKLIB lapack REQUIRED)
find_library(BLASLIB openblas REQUIRED)

message(STATUS "LAPACK FOUND: " ${LAPACKLIB})

CUDA_ADD_LIBRARY(
  gpuFunctions
#  SHARED
  gpuFunctionsImpl.cu
)

set (SHARED_OR_NOT SHARED)

set_property(TARGET gpuFunctions PROPERTY CUDA_ARCHITECTURES 52 60 70 72 75)

CUDA_ADD_LIBRARY(
  cudasirecon
#  ${SHARED_OR_NOT}
  cudaSirecon.cpp
  boostfs.cpp
)

CUDA_ADD_EXECUTABLE(cudaSireconDriver cudaSireconDriver.cpp)
CUDA_ADD_LIBRARY(pycudasirecon SHARED pyInterface.cpp)

add_dependencies(cudaSireconDriver cudasirecon gpuFunctions Buffer)
add_dependencies(pycudasirecon cudasirecon gpuFunctions Buffer)

target_link_libraries(cudasirecon gpuFunctions Buffer ${TIFF_LIBRARIES} ${LAPACKLIB} ${BLASLIB}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
)

if (IMLIB)
  add_definitions(-DMRC)
  target_link_libraries(cudasirecon ${IMLIB} ${IVELIB})
endif()

target_link_libraries(cudaSireconDriver cudasirecon)
target_link_libraries(pycudasirecon cudasirecon)

if(WIN32)
  if ( ${MSVC_VERSION} GREATER 1800)
    target_link_libraries(cudaSireconDriver ${LEG_STDIO})
  endif()
else()
  target_link_libraries(cudasirecon X11 rt)
endif(WIN32)

CUDA_ADD_CUFFT_TO_TARGET(cudasirecon)
CUDA_ADD_CUFFT_TO_TARGET(cudaSireconDriver)

# added for make install to work in conda
set(HEADERS
  cudaSirecon.h
  cudaSireconImpl.h
  gpuFunctions.h
)

install(
  TARGETS cudaSireconDriver cudasirecon
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(FILES ${HEADERS} DESTINATION include)
