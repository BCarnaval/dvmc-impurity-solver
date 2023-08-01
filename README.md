<div align="center">

# dvmc-impurity-solver

This package computes the Variational Monte Carlo ground state
and its associated Green function. This package is an extension
of the dVMC package published in
[arxiv:1912.09960](https://arxiv.org/abs/1912.09960)
and [PhysRevX.10.041023](https://doi.org/10.1103/PhysRevX.10.041023)
which is itself based on the original mVMC open source package
[source](https://github.com/issp-center-dev/mVMC)
and [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).

![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

![LICENSE](https://img.shields.io/github/license/BCarnaval/DynamicalVMC?color=blue&style=for-the-badge) ![release](https://img.shields.io/github/v/tag/BCarnaval/DynamicalVMC?color=%23FF7F50&style=for-the-badge)

</div>

It reuses most of the previous source code. Please refer to
the documentation of mVMC in the open source package
as most will not be covered here. Here, we will cover the
difference and new functions not contained in the original mVMC
and dVMC packages. This code is a direct extension of the dvmc code availabe in
[PhysRevX.10.041023](https://doi.org/10.1103/PhysRevX.10.041023).
It has essentially the same behavior but can now relax the
constraint of translation invariance and periodic boundary condition.
It also offers the possibility to include bath orbitals which
are treated differently (no calcultion of the Green function
on these orbitals).

This is very important to be able to use this Hubbard model
solution as an impurity solver that can be conjugated with
CPT, CDMFT, etc.

> #### Note
>
> The present authors (Maxime Charlebois, Peter Rosenberg and Antoine de Lagrave)
> only work on these changes and not on the original mVMC
> package. You can find the original `README` of both previous implementations
> in the `./doc` directory.

# Table of contents

- [Requirements](#requirements)

  - [Compilation](#compilation)
  - [Parallelization](#parallelization)
  - [LAPACK and BLAS](#lapack-and-blas)
  - [Python](#python)
  - [PyQCM](#pyqcm)

- [Installation](#installation)

- [Usage](#usage)

- [Authors](#authors)

- [License](#license)

# Requirements

## Compilation

You must have installed at least one if the following software to compile the C code

- [icc](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html) version >= 2021.9.0

- [gcc](https://gcc.gnu.org/) version >= 13.1.0

- [clang](https://clang.llvm.org/) version >= 13.0.0

and a Fortran compiler

- [ifort](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.3uuywf) (not tested)

- [gfortran](https://gcc.gnu.org/wiki/GFortran) version >= 13.1.0

In order to compile and link properly the code base, the software uses
CMake. Most of this project's directories
contain a specific `CMakelists.txt` with the right instructions.

- [CMake](https://cmake.org/) version >= 3.12

## Parallelization

The program is designed to be executed using supercomputers and is optimized
using parallel libraries such as [OpenMP](https://www.openmp.org/) and
OpenMPI. OpenMP is a compiler dependency, so it does not need to be separately
installed as it is included with the C/C++ compiler. However, this is not the
case for OpenMPI, it must be installed separately

- [OpenMPI](https://www.open-mpi.org/) version >= 4.1.5

## LAPACK and BLAS

LAPACK (Linear Algebra PACKage) is a software
library for numerical computation that provides routines for solving linear
algebra problems, such as linear system solving, eigenvalue and eigenvector
computations, and matrix factorizations.

- [LAPACK](https://www.netlib.org/lapack/) version >= 3.11

BLAS (Basic Linear Algebra Subprograms), on the
other hand, is a library of basic functions for linear algebra operations, such
as matrix multiplication, vector operations, and discrete Fourier transform
operations.

- [BLAS](https://www.netlib.org/blas/) version >= 0.3.23

These two libraries are often used together to achieve high-performance linear
algebra computations on modern computers.

## Python

The interface tools of this project found inside `./tool/dvmc/` are written
in Python 3. The latest **stable** version for Python 3 is 3.10.4 and is totally
compatible with [PyQCM](#pyqcm) which needs a Python version >= 3.7.

- [Python](https://www.python.org/) version >= 3.7

## PyQCM

> PyQCM is a python module that interfaces with a library written in C++ : qcm.
> This library provide a collection of functions that help implement quantum
> cluster methods. Specifically, it provides an exact diagonalization solver for
> small clusters on which a Hubbard-like model is defined and provides functions
> to define infinite-lattice models and to embed the clusters into the lattice
> via Cluster Pertrubation Theory (CPT). Methods like the Variational Cluster
> Approximation (VCA) and Cluster Dynamical Mean Field Theory (CDMFT) are then
> imple- mented from qcm by the `pyqcm` module, which is written in Python only.
>
> --David Sénéchal
>
> (<https://qcm-wed.readthedocs.io/en/latest/intro.html#what-is-pyqcm>)

- [PyQCM](https://bitbucket.org/dsenechQCM/qcm_wed/src/master/) version >= 2.2.1

# Installation

To compile and install everything, follow these instructions from root directory:

1. Making the conventional CMake working directory

```shell
mkdir build && cd build
```

2. Telling CMake to use specified configuration file based on user's compiler and OS. All
   possible configurations can be found inside the `./config/` directory

```shell
cmake ..
```

or

```shell
cmake .. -DCONFIG=<chosen_config_file.cmake>
```

to use a specifig configuration as explained above (see `./config/` directory).

3. Compiling and linking the code with specified compiler settings

```shell
make
```

4. Installing the binaries and dVMC CLI (Command Line Interface)

```shell
make install
```

It is recommended to consult the original mVMC documentation to learn
more about dependencies and parameters names and idea. Note however that
most of the cases in mVMC are not covered by the present code, as
stated above.

# Usage

## Samples

Detailed usage is not covered here. Instead many examples can be found
in the `./samples/` subdirectory. The easiest way to understand how to
use it is to run these examples. Go there to read the README.md and run the
few examples.

## dvmc CLI (Command Line Interface)

Most of the usage information is contained inside the `dvmc` command line tool
installed by CMake inside `$HOME/.local/bin` on your system. If the command

```shell
dvmc --help
```

doesn't work on your machine, you should add this directory to your `$PATH` so it can
be found as a command. To temporarily add this directory to your `$PATH` you can use the
command

```shell
export PATH="$HOME/.local/bin:$PATH"
```

By running this command, your shell will be able to access the scripts in `$HOME/.local/bin`
until the current session is ended. To add it every time a shell session is
openned, you must add the previous line to the `$HOME/.bashrc` file.

# Authors

**Maxime Charlebois**, **Peter Rosenberg** and **Antoine de Lagrave**.
Please email us at

- <maxime.charlebois@uqtr.ca>
- <peter.rosenberg@usherbrooke.ca>
- <antoine.de.lagrave@usherbrooke.ca>,

if you have any question on the code. This complete the documentation of dVMC.

# License

GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).
See `LICENSE` file for the details.
