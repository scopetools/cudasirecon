include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    "${CMAKE_CURRENT_SOURCE_DIR}/../Buffers")

find_package(LAPACK REQUIRED)

add_library(gpuFunctions gpuFunctionsImpl.cu)
add_library(cudasirecon cudaSirecon.cpp boostfs.cpp)
add_executable(cudaSireconDriver cudaSireconDriver.cpp)

cuda_add_cufft_to_target(cudasirecon)
cuda_add_cufft_to_target(cudaSireconDriver)

set_property(TARGET gpuFunctions PROPERTY MSVC_RUNTIME_LIBRARY "")
set_property(TARGET cudasirecon PROPERTY MSVC_RUNTIME_LIBRARY "")
set_property(TARGET cudaSireconDriver PROPERTY MSVC_RUNTIME_LIBRARY "")

if(BUILD_MRC)
  target_link_libraries(cudasirecon ${IMLIB} ${IVELIB})
endif()

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
  ${Boost_FILESYSTEM_LIBRARY})

target_link_libraries(cudaSireconDriver cudasirecon)

if(WIN32)
  if(${MSVC_VERSION} GREATER 1800)
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
  ARCHIVE DESTINATION lib)

install(FILES ${HEADERS} DESTINATION include)
