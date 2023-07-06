# Include guard
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(
    FATAL_ERROR
      "CMake shouldn't be executed for 'src' subdirectory, but for the top directory of UHF."
  )
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)

set(SOURCES_UHF
    UHFmain.c
    output.c
    cal_energy.c
    green.c
    makeham.c
    diag.c
    initial.c
    matrixlapack.c
    readdef.c)

include_directories(../sfmt)
add_definitions(-DMEXP=19937)

set(SOURCES_sfmt ../sfmt/SFMT.c)

add_executable(UHF ${SOURCES_UHF} ${SOURCES_sfmt})
target_link_libraries(UHF ${LAPACK_LIBRARIES} m)
