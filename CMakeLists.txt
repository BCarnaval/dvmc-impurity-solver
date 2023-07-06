# Include guard
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(mVMC NONE)

set(OPENMP_LIBRARIES "/usr/local/opt/libomp/lib")
set(OPENMP_INCLUDES "/usr/local/opt/libomp/include")

option(USE_SCALAPACK "Use Scalapack" OFF)

add_definitions(-D_mVMC)
if(CONFIG)
  message(STATUS "Loading configuration file: "
                 ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
enable_language(C Fortran)

set(MPI_HOME "/opt/openmpi")

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG} -ggdb -Wall")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran} ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  if("${CMAKE_C_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
    set(OMP_FLAG_Intel "-openmp")
  else()
    set(OMP_FLAG_Intel "-qopenmp")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OMP_FLAG_Intel}")
else()

  option(USE_OpenMP"Use OpenMP to enamble <omp.h>" ON)
  if(CMAKE_C_COMPILER_ID MATCHES "Clang")

    set(OpenMP_C_LIB_NAMES "libomp")
    set(OpenMP_CXX_LIB_NAMES "libomp")
    set(OpenMP_C_FLAGS
        "-Xclang -fopenmp -I/usr/local/opt/libomp/include/ -L /usr/local/opt/libomp/lib -lomp"
    )
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/include")
    set(OpenMP_omp_LIBRARY "/usr/local/lib/libomp.dylib")
    set(OpenMP_libomp_LIBRARY "/usr/local/lib/libomp.dylib")

  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
  endif()

  find_package(OpenMP REQUIRED)

  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif(OPENMP_FOUND)

endif()

find_package(MPI COMPONENTS C)
if(MPI_FOUND)
  message(${MPI_C_INCLUDE_PATH})
  include_directories(${MPI_C_INCLUDE_PATH})
  add_definitions(-D_mpi_use)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
else()
  message(FATAL_ERROR "MPI not found, please verify MPI installation PATH.")
endif(MPI_FOUND)

find_package(LAPACK)
if(USE_SCALAPACK MATCHES OFF)
  if(LAPACK_FOUND)
    add_definitions(-D_lapack)
  else()
    message(
      FATAL_ERROR "LAPACK not found, please verify LAPACK installation PATH.")
  endif(LAPACK_FOUND)
endif()

add_subdirectory(src/ComplexUHF)
add_subdirectory(src/pfapack)
add_subdirectory(src/mVMC)
add_subdirectory(tool)
add_subdirectory(tool/dvmc)
add_subdirectory(tool/dvmc_cli)

# Build and enable tests testing setup enable_testing() must be called in the
# top-level CMakeLists.txt before any add_subdirectory() is called.
option(Testing "Enable testing" OFF)

if(Testing)
  enable_testing(test)
  add_subdirectory(test)
endif()
