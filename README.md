

# dvmc-impurity-solver

This package computes the Variational Monte Carlo ground state
and its associated Green function for the Hubbard model. Details 
on the physics associated to this package can be found at:
[arxiv.org:2307.15738](https://arxiv.org/abs/2307.15738).

This package is an extension of the dVMC package published in
[PhysRevX.10.041023](https://doi.org/10.1103/PhysRevX.10.041023) - ([arxiv:1912.09960](https://arxiv.org/abs/1912.09960)) which is itself based on the 
original mVMC open source [mVMC package](https://dx.doi.org/10.17632/xhgyp6ncvt.1) from:
[ComputPhysCommun.08.014](https://doi.org/10.1016/j.cpc.2018.08.014) - ([arXiv:1711.11418](https://arxiv.org/abs/1711.11418))

![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white) ![Fortran](https://img.shields.io/badge/Fortran-%23734F96.svg?style=for-the-badge&logo=fortran&logoColor=white) ![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) ![Shell Script](https://img.shields.io/badge/shell_script-%23121011.svg?style=for-the-badge&logo=gnu-bash&logoColor=white)

![LICENSE](https://img.shields.io/github/license/BCarnaval/DynamicalVMC?color=blue&style=for-the-badge) ![release](https://img.shields.io/github/v/tag/BCarnaval/DynamicalVMC?color=%23FF7F50&style=for-the-badge)

</div>

This package has essentially the same behavior as the code in [PhysRevX.10.041023](https://doi.org/10.1103/PhysRevX.10.041023) but can now relax the constraint of translation invariance and periodic boundary  condition. It also offers the possibility to include bath orbitals  which are treated differently (no calcultion of the Green function on these orbitals). This new feature is very important to be able to use this Hubbard impurity solver solution as an impurity solver that can be conjugated with CPT, CDMFT, etc.

> #### Note
>
> The present authors (Maxime Charlebois, Peter Rosenberg and Antoine de Lagrave)
> only work on these changes and not on the original mVMC package. 
> You can find the original `README` of both original mVMC and dVMC 
> packages in the `./doc` directory.

# Table of contents

- [References](#references)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)
- [License](#license)

# References

You can look for the latest developments of this package at [github.com/BCarnaval/dvmc-impurity-solver](https://github.com/BCarnaval/dvmc-impurity-solver).

You are free to use this code as long as you respect the License terms. If you use it or learn from it, we ask that you cite the following article:
[arxiv.org:2307.15738](https://arxiv.org/abs/2307.15738): CDMFT-dVMC article, to be published in PhysRevB,

and potentially these articles:
[PhysRevB.106.245132](https://doi.org/10.1103/PhysRevB.106.245132) - ([arxiv:2209.08092](https://arxiv.org/abs/2209.08092)) : CPT-dVMC article.
[PhysRevX.10.041023](https://doi.org/10.1103/PhysRevX.10.041023) - ([arxiv:1912.09960](https://arxiv.org/abs/1912.09960)) : original dVMC article and software library.
[ComputPhysCommun.08.014](https://doi.org/10.1016/j.cpc.2018.08.014) - ([arXiv:1711.11418](https://arxiv.org/abs/1711.11418)): original VMC article and software library.

# Dependencies

In summary it requires C and Fortran compilers, cmake, openmpi, lapack-blas, python3 and pyqcm. Most of dependencies can be installed using these simple commands:

### Linux (Ubuntu)
```shell
sudo apt update
sudo apt install git build-essential gcc gfortran make cmake libblas-dev liblapack-dev
sudo apt install openmpi-bin openmpi-doc libopenmpi-dev checkinstall
```

### MacOS
```shell  
xcode-select --install
brew install git cmake open-mpi openblas lapack
```
After that, a compatible version of PyQCM must be installed according to the procedure detailed at the end of `./doc/INSTALL_DEPENDENCIES.md`. More information about the **installation procedures** and **compatible versions** of dependencies are detailed in this file.

# Installation

1. Clone the repository [dvmc-impurity-solver](https://github.com/BCarnaval/dvmc-impurity-solver) on your machine

```shell
git clone https://github.com/BCarnaval/dvmc-impurity-solver
```

2. Making the conventional CMake working directory

```shell
cd dvmc-impurity-solver
mkdir build && cd build
```

3. Telling CMake to use specified configuration file based on user's compiler and OS. All
   possible configurations can be found inside the `./config/` directory

```shell
cmake ..
```

or

```shell
cmake .. -DCONFIG=<chosen_config_file.cmake>
```

to use a specifig configuration as explained above (see `./config/` directory).

4. Compiling and linking the code with specified compiler settings and installing the binaries and dVMC CLI 

```shell
make install
```

5. Modifying `$PYTHONPATH` variable to make `dvmc.py` accessible as a module

```shell
export PYTHONPATH="$HOME/.local/share/dvmc:$PYTHONPATH"
```

(Note: Add this last line to your `.bashrc` to make dVMC always accessible)

It is recommended to consult the original mVMC documentation to learn
more about dependencies and parameters names. Note however that
most of the cases present in the originale mVMC documentation are 
not covered by the present code.

# Usage

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

## Samples

Detailed usage is not covered here. Instead many examples can be found
in the `./samples/` subdirectory. The easiest way to understand how to
use it is to run these examples. Go there to read the `README.md` of the 
examples that increase in complexity from `01` to `05`.

# Authors

All changes from the original [mVMC package](https://dx.doi.org/10.17632/xhgyp6ncvt.1) were done by
**Maxime Charlebois**, **Peter Rosenberg** and **Antoine de Lagrave**.
Please email us if you have any question on this version of the code.
- <maxime.charlebois@uqtr.ca>
- <peter.rosenberg@usherbrooke.ca>
- <antoine.de.lagrave@usherbrooke.ca>,

# License

GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)).
See `LICENSE` file for the details.
