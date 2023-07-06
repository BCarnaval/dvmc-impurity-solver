# for Intel Compiler

set(CMAKE_C_COMPILER
    "icc"
    CACHE STRING "" FORCE)

set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wformat -Werror=format-security")

set(CMAKE_C_FLAGS_RELEASE
    "-Wno-unknown-pragmas -O3 -DNDEBUG -DHAVE_SSE2"
    CACHE STRING "" FORCE)

set(CMAKE_Fortran_COMPILER
    "gfortran"
    CACHE STRING "" FORCE)

set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -DNDEBUG -DHAVE_SSE2"
    CACHE STRING "" FORCE)

# for Intel MKL
set(BLA_VENDOR OpenBLAS)

list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/lapack")
list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/openblas")
