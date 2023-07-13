# Include guard
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(
    FATAL_ERROR
      "CMake shouldn't be executed for 'src' subdirectory, but for the top directory of mVMC."
  )
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
add_definitions(-DMEXP=19937)

set(SOURCES_vmcmain vmcmain.c splitloop.c)

set(SOURCES_sfmt ../sfmt/SFMT.c)

add_library(
  StdFace STATIC
  ../StdFace/ChainLattice.c
  ../StdFace/HoneycombLattice.c
  ../StdFace/SquareLattice.c
  ../StdFace/StdFace_main.c
  ../StdFace/StdFace_ModelUtil.c
  ../StdFace/TriangularLattice.c
  ../StdFace/Ladder.c
  ../StdFace/Kagome.c
  ../StdFace/Orthorhombic.c
  ../StdFace/Pyrochlore.c
  ../StdFace/Wannier90.c
  ../StdFace/FCOrtho.c)

link_directories("$ENV{HOME}/.local/share/dvmc/lib")

add_executable(dvmcdry.out vmcdry.c)
target_link_libraries(dvmcdry.out StdFace m)
add_executable(dvmc.out ${SOURCES_vmcmain} ${SOURCES_sfmt})
if(MPI_FOUND)
  target_link_libraries(dvmc.out pfapack StdFace ${LAPACK_LIBRARIES} m
                        ${MPI_C_LIBRARIES})
else()
  target_link_libraries(dvmc.out pfapack StdFace ${LAPACK_LIBRARIES} m)
endif(MPI_FOUND)

if(USE_SCALAPACK)
  string(REGEX REPLACE "-L[ ]+" "-L" sc_libs ${SCALAPACK_LIBRARIES})
  string(REGEX REPLACE "[ ]+" ";" sc_libs ${sc_libs})
  foreach(sc_lib IN LISTS sc_libs)
    target_link_libraries(dvmc.out ${sc_lib})
  endforeach(sc_lib)
  message(STATUS "SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
endif(USE_SCALAPACK)

install(TARGETS dvmcdry.out RUNTIME DESTINATION $ENV{HOME}/.local/share/dvmc)
install(TARGETS dvmc.out RUNTIME DESTINATION $ENV{HOME}/.local/share/dvmc)
add_definitions(-D_mVMC)
