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
    include_directories( ${TIFF_INCLUDE_DIR} )
    link_directories ( ${TIFF_LIBRARY} )
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
  find_package(LAPACK)
  find_library(IMLIB libimlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
  find_library(IVELIB libive PATHS ${PRIISM_LIB_PATH} REQUIRED)

  # if ( ${MSVC_VERSION} GREATER 1800)
  #   find_library(LEG_STDIO legacy_stdio_definitions)
  #   message(STATUS "adding legacy stdio definitions for >VS14.0 " ${LEG_STDIO})
  # endif()
  message(STATUS "LAPACK: " ${LAPACK_LIBRARIES})

else()
  find_library(IMLIB imlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
  find_library(IVELIB ive PATHS ${PRIISM_LIB_PATH} REQUIRED)
endif()

message(STATUS "IMLIB: " ${IMLIB})
message(STATUS "IVELIB: " ${IVELIB})
message(STATUS "LEG_STDIO: " ${LEG_STDIO})


find_package(CUDA)
message(STATUS "CUDA: " ${CUDA_CUFFT_LIBRARIES})


CUDA_ADD_LIBRARY(
  gpuFunctions
#  SHARED
  gpuFunctionsImpl.cu
  )

if (NOT APPLE)
set (SHARED_OR_NOT SHARED)
endif()

CUDA_ADD_LIBRARY(
  cudaSirecon
#  ${SHARED_OR_NOT}
  cudaSirecon.cpp
  boostfs.cpp
  )

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
  target_link_libraries(
    cudaSirecon
    gpuFunctions
    Buffer
    ${IMLIB}
    ${IVELIB}
    ${LAPACK_LIBRARIES}
    ${TIFF_LIBRARIES}
    )

  if ( ${MSVC_VERSION} GREATER 1800)
    target_link_libraries(
      cudaSireconDriver
      cudaSirecon
      gpuFunctions
      Buffer
      ${LEG_STDIO}
      ${IMLIB}
      ${IVELIB}
      ${LAPACK_LIBRARIES}
      )
  else()
	target_link_libraries(
      cudaSireconDriver
      cudaSirecon
      gpuFunctions
      Buffer
      ${IMLIB}
      ${IVELIB}
      ${LAPACK_LIBRARIES}
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
