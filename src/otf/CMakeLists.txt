if(WIN32)
  set(PLATFORM win64)
else()
  set(PLATFORM linux64)
endif(WIN32)

include_directories("${CMAKE_SOURCE_DIR}/cudaSirecon")  # for CImg.h
add_executable(makeotf radialft.cpp)

if (BUILD_MRC)
  set(PRIISM_LIB_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/LIB")
  set(PRIISM_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/IVE/${PLATFORM}/INCLUDE")
  include_directories(${PRIISM_INCLUDE_PATH})
  link_directories(${PRIISM_LIB_PATH})

  find_library(IMLIB NAMES imlib libimlib PATHS ${PRIISM_LIB_PATH} REQUIRED)
  find_library(IVELIB NAMES ive libive PATHS ${PRIISM_LIB_PATH} REQUIRED)
  add_definitions(-DMRC)
  target_link_libraries(makeotf ${IMLIB} ${IVELIB})
endif()



if (WIN32)
  find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB)
  # include_directories( ${FFTW_INCLUDE_DIRS} )
  if (${MSVC_VERSION} GREATER 1800)
    target_link_libraries(makeotf ${LEG_STDIO})
  endif()
else()
  find_package(FFTW REQUIRED COMPONENTS FLOAT_THREADS_LIB)
endif()
  
message(STATUS "FFTW3 libraries: " ${FFTW_LIBRARIES})
target_link_libraries(makeotf ${FFTW_LIBRARIES} ${TIFF_LIBRARIES})


install(TARGETS makeotf RUNTIME DESTINATION bin)


if (BUILD_OTF_VIEWER)
  add_executable(otfviewer otfviewer.cpp)
  target_link_libraries(otfviewer ${IMLIB} ${IVELIB} ${TIFF_LIBRARIES})
  
  if (WIN32)
    if ( ${MSVC_VERSION} GREATER 1800)
      target_link_libraries(otfviewer ${LEG_STDIO})
    endif()
  else()
    target_link_libraries(otfviewer X11)
  endif()

  install(TARGETS otfviewer RUNTIME DESTINATION bin)
endif()
