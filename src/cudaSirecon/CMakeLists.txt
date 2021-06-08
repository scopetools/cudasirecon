if(WIN32)
  link_directories ( ${PRIISM_LIB_PATH} )
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/win64/INCLUDE")
else()
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/linux64/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/linux64/INCLUDE")
endif(WIN32)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers"
  ${PRIISM_INCLUDE_PATH}
)

find_library(IMLIB NAMES imlib libimlib PATHS ${PRIISM_LIB_PATH})
find_library(IVELIB NAMES ive libive PATHS ${PRIISM_LIB_PATH})
find_package(LAPACK REQUIRED)

message(STATUS "IMLIB: " ${IMLIB})
message(STATUS "IVELIB: " ${IVELIB})
message(STATUS "LEG_STDIO: " ${LEG_STDIO})
message(STATUS "LAPACK_FOUND: " ${LAPACK_FOUND})

CUDA_ADD_LIBRARY(gpuFunctions gpuFunctionsImpl.cu)
set_property(TARGET gpuFunctions PROPERTY CUDA_ARCHITECTURES 52 60 70 72 75)

CUDA_ADD_LIBRARY(cudasirecon cudaSirecon.cpp boostfs.cpp)
CUDA_ADD_EXECUTABLE(cudaSireconDriver cudaSireconDriver.cpp)
set_property(TARGET cudaSireconDriver PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")

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

if (IMLIB)
  add_definitions(-DMRC)
  target_link_libraries(cudasirecon ${IMLIB} ${IVELIB})
endif()

target_link_libraries(cudaSireconDriver cudasirecon)

if(WIN32)
  if ( ${MSVC_VERSION} GREATER 1800)
    target_link_libraries(cudaSireconDriver ${LEG_STDIO})
  endif()
else()
  target_link_libraries(cudasirecon rt)
endif(WIN32)

CUDA_ADD_CUFFT_TO_TARGET(cudasirecon)
CUDA_ADD_CUFFT_TO_TARGET(cudaSireconDriver)

# added for make install to work in conda
set(HEADERS cudaSirecon.h cudaSireconImpl.h gpuFunctions.h)

install(
  TARGETS cudaSireconDriver cudasirecon
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(FILES ${HEADERS} DESTINATION include)
