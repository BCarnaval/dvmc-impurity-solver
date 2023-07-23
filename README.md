<div align="center">

# DynamicalVMC

This code is based on the original mVMC open source package
[source](https://github.com/issp-center-dev/mVMC)
and [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).
It reuse most of the previous source code. Please refer to
the documentation of mVMC in the open source package
as most will not be covered here. Here, we will cover the
difference and new functions not contained in the original mVMC.

This package is not as complete as the mVMC package.
It is only a prototype used to calculate the "dvmc" method implemented
and published in [arxiv:1912.09960](https://arxiv.org/abs/1912.09960).
However, it is very useful and easy to reproduce some date in this
publication.

![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

![LICENSE](https://img.shields.io/github/license/BCarnaval/DynamicalVMC?color=blue&style=for-the-badge) ![release](https://img.shields.io/github/v/tag/BCarnaval/DynamicalVMC?color=%23FF7F50&style=for-the-badge)

</div>

This code implements a new feature where you can calculate the
frequency dependent Green function (option NVMCCalMode = 3).
For now there is a number of condition on this
new calculation mode. It is restricted to:

- 1D or 2D model (easy to generalize to 3D, to be done)
- real number calculation
- OrbitalAntiParallel (not OrbitalGeneral or OrbitalParallel)
- Hubbard model

The present authors (Maxime Charlebois, Peter Rosenberg and Antoine de Lagrave) only work on these
changes and not on the original mVMC package. You can find the original README
and authors of mVMC at the end of this document.

# Table of contents

- [Requirements](#requirements)

  - [Compilation](#compilation)
  - [Parallelization](#parallelization)
  - [LAPACK and BLAS](#lapack-and-blas)
  - [Python](#python)

- [Installation](#installation)

- [Usage](#usage)

- [Details](#details)

- [Original mVMC README](#original-mvmc-readme)

# Requirements

## Compilation

The project being written in C, a C compiler (`icc`, `gcc`, `clang`) is then
required to use this software. It also has Fortran dependencies:
[pfapack](https://michaelwimmer.org/downloads.html) which also need to be
compiled using Fortran compiler (`ifort`, `gfortran`).

In order to compile and link properly the code base, the software uses
[CMake](https://cmake.org/) version >= 3.12. Most of this project's directories
contain a specific `CMakelists.txt` with the right instructions. This tool can be
installed via

### MacOS

```shell
brew install cmake
```

### Linux (Ubuntu)

```shell
sudo apt install cmake
```

## Parallelization

The program is designed to be executed using supercomputers and is optimized
using parallel libraries such as [OpenMP](https://www.openmp.org/) and
[OpenMPI](https://www.open-mpi.org/). OpenMP is a compiler dependency, so it
does not need to be separately installed as it is included with the C/C++
compiler. However, this is not the case for OpenMPI, and it can be installed
and verified using the following commands:

### MacOS

```shell
brew install open-mpi
```

### Linux (Ubuntu)

```shell
sudo apt install openmpi-bin openmpi-doc libopenmpi-dev
```

## LAPACK and BLAS

[LAPACK](https://www.netlib.org/lapack/) (Linear Algebra PACKage) is a software
library for numerical computation that provides routines for solving linear
algebra problems, such as linear system solving, eigenvalue and eigenvector
computations, and matrix factorizations.

[BLAS](https://www.netlib.org/blas/) (Basic Linear Algebra Subprograms), on the
other hand, is a library of basic functions for linear algebra operations, such
as matrix multiplication, vector operations, and discrete Fourier transform
operations.

These two libraries are often used together to achieve high-performance linear
algebra computations on modern computers.

### MacOS

```shell
brew install openblas lapack
```

### Linux (Ubuntu)

```shell
sudo apt install libblas-dev liblapack-dev
```

## Python

The interface tools of this project found inside `./tool/dvmc/` are written
in Python 3. The latest **stable** version for Python 3 is 3.10.4 and is totally
compatible with [PyQCM](https://bitbucket.org/dsenechQCM/qcm_wed/src/master/).
The installation could be done by

### MacOS

```shell
brew install python@3.10
```

### Linux (Ubuntu)

```shell
sudo apt install wget libncursesw5-dev libssl-dev libsqlite3-dev tk-dev \
libgdbm-dev libc6-dev libbz2-dev libffi-dev zlib1g-dev
...
sudo apt install python3.10 python3-pip
```

# Installation

To compile and install everything, follow these instructions from root directory:

1. Making the conventional CMake working directory

```shell
mkdir build && cd build
```

2. Telling CMake to use specified configuration file based on user's compiler and OS. All
   possible configurations can be found inside the `./config/` directory

```shell
cmake .. -DCONFIG=<chosen_config_file.cmake>
```

3. Compiling and linking the code with specified compiler settings

```shell
make
```

4. Installing the binaries and dVMC CLI symlink (Command Line Interface)

```shell
make install
```

It is recommended to consult the original mVMC documentation to learn
more about dependencies and parameters names and idea. Note however that
not all the case in mVMC are covered by the present code, as stated above.

# Usage

## Samples

Detailed usage is not covered here. Instead many examples can be found
in the `./samples/` subdirectory. The easiest way to understand how to
use it is to run these examples. Go there to read the README and run the
few examples.

## dvmc CLI

Most of the usage information is contained inside the `dvmc` command line tool
installed by CMake inside `$HOME/.local/bin` on your system. If the command

```shell
dvmc --help
```

doesn't work on your machine, you should add this directory to your `$PATH` so it can
be found as a command.

# Details

Please email us at <maxime.charlebois@uqtr.ca> or <antoine.de.lagrave@usherbrooke.ca>
if you have any question on the code. This complete the documentation of dVMC.

# Original mVMC README

A numerical solver package for a wide range of quantum lattice models based on many-variable Variational Monte Carlo method

## What is mVMC ?

mVMC (many-variable Variational Monte Carlo method)
is a software for performing the highly-accurate
variational Monte Carlo calculations
with the simple and flexible user interface.
mVMC also supports the large-scale parallelization.
For the conventional models in strongly correlated electron systems such as the Hubbard model, the Heisenberg model, and the Kondo-lattice model,
users can perform the calculation by preparing the one input files whose length is shorter than ten lines.

By using the same input file, users can perform the exact diagonalization through [HPhi](https://github.com/QLMS/HPhi/releases).
Thus, it is easy to examine the accuracy of the variational calculation for small system sizes
and to perform the calculations
for large system sizes that can not be treated
by the exact diagonalization.
A broad spectrum of users including experimental scientists is cordially welcome.

## Methods

many-variable variational Monte Carlo method

## Target models

Hubbard model, Heisenberg model, Kondo lattice model, multi-orbital Hubbard model

## Available physical quantities

specific heat, susceptibility, ground state energy, structure factors

## Requirement

- C compiler (intel, Fujitsu, GNU, etc. )
- ScaLAPACK library (intel MKL, Fujitsu, ATLAS, etc.)
- MPI library

## Install

You can install mVMC and also get a manual for mVMC from a [release note](https://github.com/issp-center-dev/mVMC/releases).

## Licence

GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).

The mVMC package is developed based on the [mVMC-mini](https://github.com/fiber-miniapp/mVMC-mini) program. The license of mVMC-mini is "The BSD 3-Clause License".

We would appreciate if you cite the following article in your research with mVMC:  
mVMC - Open-source software for many-variable variational Monte Carlo method, Takahiro Misawa, Satoshi Morita, Kazuyoshi Yoshimi, Mitsuaki Kawamura, Yuichi Motoyama, Kota Ido, Takahiro Ohgoe, Masatoshi Imada, Takeo Kato, [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).

## Tutorials

Lecture notes and sample scripts used in Hands-on
are available at [mVMC-tutorial](https://github.com/issp-center-dev/mVMC-tutorial)

## Authors

Takahiro Misawa, Satoshi Morita, Takahiro Ohgoe, Kota Ido, Yuichi Motoyama, Mitsuaki Kawamura, Kazuyoshi Yoshimi, Takeo Kato, Masatoshi Imada.
