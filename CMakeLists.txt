# -------------
# Include guard
# -------------
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(mVMC NONE)
option(USE_SCALAPACK "Use Scalapack" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
add_definitions(-D_mVMC)

# ----------------------------------------------
# Set specific configuration using -DCONFIG flag
# ----------------------------------------------
if(CONFIG)
  message(STATUS "Loading configuration file: "
                 ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)
enable_language(C Fortran)

# ---------------------------------------------------------
# Set build type (Debug or Release) 'Release being default'
# ---------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" MATCHES "Release")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Type of build" FORCE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran} ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(CMAKE_BUILD_TYPE
      "Debug"
      CACHE STRING "Type of build" FORCE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG} -ggdb -Wall")
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

# -------------------------------
# Compiler specific configuration
# -------------------------------
if(CMAKE_C_COMPILER STREQUAL "icc")
  set(OpenMP_C_LIB_NAMES "libomp" "omp")
  set(OpenMP_omp_LIBRARY "/usr/local/lib/libomp.dylib")
  set(OpenMP_libomp_LIBRARY "/usr/local/lib/libomp.dylib")
  set(OpenMP_C_FLAGS
      "-I/usr/local/opt/libomp/include/ -L /usr/local/opt/libomp/lib")
  set(OMP_FLAG_Intel "-qopenmp")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OMP_FLAG_Intel}")
elseif(CMAKE_C_COMPILER MATCHES "clang")
  set(OpenMP_C_LIB_NAMES "libomp" "omp")
  set(OpenMP_C_FLAGS
      "-Xclang -fopenmp -I/usr/local/opt/libomp/include/ -L /usr/local/opt/libomp/lib -lomp"
  )
  set(OpenMP_omp_LIBRARY "/usr/local/lib/libomp.dylib")
  set(OpenMP_libomp_LIBRARY "/usr/local/lib/libomp.dylib")
endif()

# ------------------------------------------
# OpenMP package search (fails if not found)
# ------------------------------------------
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# -------------------------------------------
# OpenMPI package search (fails if not found)
# -------------------------------------------
find_package(MPI COMPONENTS C)
if(MPI_FOUND)
  message(${MPI_C_INCLUDE_PATH})
  include_directories(${MPI_C_INCLUDE_PATH})
  add_definitions(-D_mpi_use)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
else()
  message(FATAL_ERROR "MPI not found, please verify MPI installation PATH.")
endif(MPI_FOUND)

# ------------------------------------------
# LAPACK package search (fails if not found)
# ------------------------------------------
find_package(LAPACK)
if(USE_SCALAPACK MATCHES OFF)
  if(LAPACK_FOUND)
    add_definitions(-D_lapack)
  else()
    message(
      FATAL_ERROR "LAPACK not found, please verify LAPACK installation PATH.")
  endif(LAPACK_FOUND)
endif()

# ----------------------------------------
# Adding sub-CMakeLists.txt to be compiled
# ----------------------------------------
add_subdirectory(src/ComplexUHF)
add_subdirectory(src/pfapack)
add_subdirectory(src/mVMC)
add_subdirectory(tool)
add_subdirectory(tool/dvmc)
add_subdirectory(tool/dvmc_cli)

# ---------------------------------------------------------------------------
# Build and enable tests testing setup enable_testing() must be called in the
# top-level CMakeLists.txt before any add_subdirectory() is called.
# ---------------------------------------------------------------------------
option(Testing "Enable testing" OFF)

if(Testing)
  enable_testing(test)
  add_subdirectory(test)
endif()
