include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers"
)

if(WIN32)
  set(PLATFORM win64)
else()
  set(PLATFORM linux64)
endif(WIN32)

find_package(LAPACK REQUIRED)

CUDA_ADD_LIBRARY(gpuFunctions gpuFunctionsImpl.cu)
CUDA_ADD_LIBRARY(cudasirecon cudaSirecon.cpp boostfs.cpp)
CUDA_ADD_EXECUTABLE(cudaSireconDriver cudaSireconDriver.cpp)

CUDA_ADD_CUFFT_TO_TARGET(cudasirecon)
CUDA_ADD_CUFFT_TO_TARGET(cudaSireconDriver)

set_property(TARGET gpuFunctions PROPERTY MSVC_RUNTIME_LIBRARY "")
set_property(TARGET cudasirecon PROPERTY MSVC_RUNTIME_LIBRARY "")
set_property(TARGET cudaSireconDriver PROPERTY MSVC_RUNTIME_LIBRARY "")

if (BUILD_MRC)
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/INCLUDE")
  include_directories(${PRIISM_INCLUDE_PATH})
  link_directories(${PRIISM_LIB_PATH})

  find_library(IMLIB NAMES imlib libimlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
  find_library(IVELIB NAMES ive libive PATHS ${PRIISM_LIB_PATH} REQUIRED)
  add_definitions(-DMRC)
  target_link_libraries(cudasirecon ${IMLIB} ${IVELIB})
endif ()

message(STATUS "IMLIB: " ${IMLIB})
message(STATUS "IVELIB: " ${IVELIB})
message(STATUS "LEG_STDIO: " ${LEG_STDIO})
message(STATUS "LAPACK_FOUND: " ${LAPACK_FOUND})

add_dependencies(cudaSireconDriver cudasirecon gpuFunctions Buffer)
target_link_libraries(
  cudasirecon
  gpuFunctions
  Buffer
  ${TIFF_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BOOST_LIBRARIES}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
)

target_link_libraries(cudaSireconDriver cudasirecon)

if(WIN32)
  if ( ${MSVC_VERSION} GREATER 1800)
    target_link_libraries(cudaSireconDriver ${LEG_STDIO})
  endif()
else()
  target_link_libraries(cudasirecon rt)
endif(WIN32)

# added for make install to work in conda
set(HEADERS cudaSirecon.h cudaSireconImpl.h gpuFunctions.h)

install(
  TARGETS cudaSireconDriver cudasirecon
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(FILES ${HEADERS} DESTINATION include)
