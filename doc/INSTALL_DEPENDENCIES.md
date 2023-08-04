# Detailed dependencies documentation

# Table of Contents

- [Package Manager](#package-manager)
- [Git](#git)
- [Compilers](#compilers)
- [Parallelization](#parallelization)
- [LAPACK and BLAS](#lapack-and-blas)
- [PyQCM](#pyqcm)

# Package Manager

To install everything required for the proper functioning of the program, it is advisable to have a package manager. This greatly facilitates the installation of dependencies and helps manage their updates. Package managers are used to easily install programs and handle their updates.

### Linux (Ubuntu)

Ubuntu already comes with a fairly convenient package manager called `apt` (Advanced Package Tool). It should already be available on the system. It is also essential to update the list of available packages so that the latest version is used during installation. This can be achieved with the following command:

```shell
apt update
```

However, the `sudo` command, which allows the execution of commands with administrative privileges, may not be present by default. You can install it using:

```shell
apt install sudo
```

Once `sudo` is installed, you can update all the installed packages on the system using the command:

```shell
sudo apt update
```

### MacOS

The easiest way to install these tools on MacOS is through [Homebrew](https://brew.sh/), the best package manager for MacOS. The following command will install Homebrew:

```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

If Homebrew is not already installed, it will prompt you during the installation process to install the `Xcode Command Line Tools`, which will be done automatically. For users who already have Homebrew installed, this step will be skipped.

> **Apple Silicon Machines**
>
> To complete the installation of Homebrew, you need to add it to the `$PATH` using the following commands:

```shell
echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
eval "$(/opt/homebrew/bin/brew shellenv)"
```

In any case, it is good practice to check the installation of Homebrew:

```shell
brew doctor
```

---

**NOTE:**

_From here onwards, we assume that the Linux system on which the user is installing the program is up to date, along with its libraries, and has the `sudo` command available._

---

# Git

To access the program and obtain the project's history, the user needs [git](https://git-scm.com/). This program allows, among other things, to track changes made to the project and collaborate with other `git` users. It is an essential tool for software development.

### Linux (Ubuntu)

```shell
sudo apt install git
```

### MacOS

```shell
brew install git
```

# Compilers

### Linux (Ubuntu)

For Linux, an equivalent to the `Xcode Command Line Tools` is a collection of compilation tools called `Ubuntu Build Essentials`. It includes:

- g++ (C++ compiler)
- gcc (C compiler)
- make (guide for compilation)
- libc6-dev (contains the 'header files' and some C libraries)
- dpkg-dev (allows packaging of applications/scripts)

To install all these tools, Ubuntu recommends the following command:

```shell
sudo apt install build-essential
```

However, these tools do not include a Fortran compiler, which is necessary for compiling the project. Therefore, you can install it using:

```shell
sudo apt install gfortran
```

Next, it is important to test the installation of these compilers. This can be done using the following commands:

```shell
g++ -v && gcc -v && make -v
```

The last tool to install, and not the least, is [CMake](https://cmake.org/). CMake allows the generation of build files like `makefile` in a general way on different operating systems, programming languages, and project types. In our case, this tool facilitates the overall compilation of the program and the integration of parallelization libraries. To install CMake on Ubuntu:

```shell
sudo apt install cmake
```

### MacOS

For working on C/C++ projects, it is strongly recommended to install the `Xcode Command Line Tools` (a collection of tools for MacOS developers such as c++, cc, clang, clang++, cpp, g++, gcc, git, ld, llvm-, make, python3.8, pip3, etc.). Install the `Xcode Command Line Tools` using:

```shell
xcode-select --install
```

or update these tools (and others recommended by the system) using:

```shell
softwareupdate -i -r
```

The last tool to install, and not the least, is [CMake](https://cmake.org/). CMake allows the generation of build files like `makefile` in a general way on different operating systems, programming languages, and project types. In our case, this tool facilitates the overall compilation of the program and the integration of parallelization libraries. To install CMake on MacOS, simply use the Homebrew package manager:

```shell
brew install cmake
```

# Parallelization

The program is designed to be executed using supercomputers, and it is optimized using parallel libraries such as [OpenMP](https://www.openmp.org/) and [OpenMPI](https://www.open-mpi.org/). OpenMP is a compiler dependency, so it does not need to be installed separately as it is included with the C/C++ compiler. However, this is not the case for OpenMPI, and it is possible to install and verify its installation using the following commands:

### Linux (Ubuntu)

```shell
sudo apt install openmpi-bin openmpi-doc libopenmpi-dev checkinstall
```

### MacOS

```shell
brew install open-mpi
```

# LAPACK and BLAS

[LAPACK](https://www.netlib.org/lapack/) (Linear Algebra PACKage) is a software library for numerical computation that provides routines for solving linear algebra problems, such as linear system solving, eigenvalue and eigenvector computations, and matrix factorizations.

[BLAS](https://www.netlib.org/blas/) (Basic Linear Algebra Subprograms), on the other hand, is a library of basic functions for linear algebra operations, such as matrix multiplication, vector operations, and discrete Fourier transform operations.

These two libraries are often used together to provide high-performance linear algebra computations on modern computers.

### Linux (Ubuntu)

Installation is done as follows:

```shell
sudo apt install libblas-dev liblapack-dev
```

### MacOS

These two libraries can be installed using Homebrew:

```shell
brew install openblas lapack
```

# PyQCM

> PyQCM is a python module that interfaces with a library written in C++: qcm. This library provides a collection of functions that help implement quantum cluster methods. Specifically, it provides an exact diagonalization solver for small clusters on which a Hubbard-like model is defined and provides functions to define infinite-lattice models and to embed the clusters into the lattice via Cluster Perturbation Theory (CPT). Methods like the Variational Cluster Approximation (VCA) and Cluster Dynamical Mean Field Theory (CDMFT) are then implemented from qcm by the `pyqcm` module, which is written in Python only.
>
> \--David Sénéchal
>
> ([https://qcm-wed.readthedocs.io/en/latest/intro.html#what-is-pyqcm](https://qcm-wed.readthedocs.io/en/latest/intro.html#what-is-pyqcm))

### Default `pip` Installation

To proceed with this type of installation, you need to have Python (>=3.7) installed. A recommended version at the time of writing this text is Python3.10 since Python3.11 is not yet stable. To install Python 3.10 on MacOS, use Homebrew:

```shell
brew install python@3.10
```

For Linux (Ubuntu), you need to install some necessary libraries for Python to function correctly:

```shell
sudo apt install wget libncursesw5-dev libssl-dev libsqlite3-dev tk-dev \
libgdbm-dev libc6-dev libbz2-dev libffi-dev zlib1g-dev
```

Then install Python and its package manager `pip` (which is not included by default in the `apt` installation):

```shell
sudo apt install python3.10 python3-pip
```

Assuming that Python and `pip` are correctly installed on the system, the PyQCM library depends on a few Python libraries that are essential to install:

```shell
pip3 install numpy scipy matplotlib
```

The user can now install PyQCM by first cloning the project to the system using `git`:

```shell
git clone https://bitbucket.org/dsenechQCM/qcm_wed.git
```

Then go to the project directory and proceed with the installation:

```shell
cd ./qcm_wed && pip install .
```

---

**NOTE:**

The user can also choose to use optional libraries during the installation to improve the performance of PyQCM. These libraries are: [CUBA](https://feynarts.de/cuba/), [PRIMME](https://www.cs.wm.edu/~andreas/software/), and [Eigen](https://gitlab.com/libeigen/eigen). To use these dependencies, provide the following directives to the `cmake` build tool during installation:

```shell
export CMAKE_ARGS="-DDOWNLOAD_CUBA=1 -DEIGEN_HAMILTONIAN=1 -DWITH_PRIMME=1 -DDOWNLOAD_PRIMME=1"
```

just before proceeding with the `pip` installation. Don't forget to install the Eigen library separately on MacOS and Linux (Ubuntu):

```shell
brew install eigen
```

```shell
sudo apt install libeigen3-dev
```

---

### Performance '_from_source_' Installation

This type of installation is more _local_, meaning that the library will not be installed within a Python environment but rather compiled directly in the source code folder. To do this, access the source code using `git`:

```shell
git clone https://bitbucket.org/dsenechQCM/qcm_wed.git
```

Next, you need to compile the library. Conventionally, with `CMake`, you proceed in a `build` directory:

```shell
cd qcm_wed && mkdir build && cd build
```

In this directory, you configure the compilation options to use high-performance external libraries: [CUBA](https://feynarts.de/cuba/), [PRIMME](https://www.cs.wm.edu/~andreas/software/), and [Eigen](https://gitlab.com/libeigen/eigen)

```shell
cmake .. -DDOWNLOAD_CUBA=1 -DEIGEN_HAMILTONIAN=1 -DWITH_PRIMME=1 -DDOWNLOAD_PRIMME=1
...
cmake --build .
...
cp ./qcm* ../pyqcm/.
```

You should manually add the library to the general Python library path of your system:

```shell
export PYTHONPATH="$HOME/qcm_wed:$PYTHONPATH"
```

This way, Python scripts executed with the default Python interpreter (excluding virtual environments) will be able to access and recognize the PyQCM modules.
