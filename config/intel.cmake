# for Intel Compiler

set(CMAKE_C_COMPILER
    "icc"
    CACHE STRING "" FORCE)

# icc compiler will be deprecated in the second half of 2023, the warning is
# disabled via: -diag-disable=10441 compilation flag.
set(CMAKE_C_FLAGS_DEBUG
    "-O0 -g -Wall -Wformat -Werror=format-security -diag-disable=10441")
set(CMAKE_C_FLAGS_RELEASE
    "-Wno-unknown-pragmas -O3 -DNDEBUG -DHAVE_SSE2 -diag-disable=10441"
    CACHE STRING "" FORCE)

set(CMAKE_Fortran_COMPILER
    "ifort"
    CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE
    "-O3 -DNDEBUG -DHAVE_SSE2"
    CACHE STRING "" FORCE)

# for Intel MKL
set(BLA_VENDOR
    "Intel10_64lp"
    CACHE STRING "" FORCE)
