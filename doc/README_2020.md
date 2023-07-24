# dVMC

Author: Maxime Charlebois

This code is based on the original mVMC open source package 
[source](https://github.com/issp-center-dev/mVMC) 
and [arXiv:1711.11418](https://arxiv.org/abs/1711.11418).
It reuse most of the previous source code. Please refer to 
the documentation of mVMC in the open source package
as most will not be covered here. Herew, we will only cover the 
difference and new functions not contained in the original mVMC.

This package is not as complete as the mVMC package.
It is only a prototype used to calculate the "dvmc" method implemented
and published in [arxiv:1912.09960](https://arxiv.org/abs/1912.09960).
However, it is very useful and easy to reproduce some date in this
publication.

This code implements a new feature where you can calculate the 
frequency dependent Green function (option NVMCCalMode = 3). 
For now there is a number of condition on this 
new calculation mode. It is restricted to:
- 1D or 2D model (easy to generalize to 3D, to be done)
- real number calculation
- OrbitalAntiParallel (not OrbitalGeneral or OrbitalParallel)
- Hubbard model

The present author (Maxime Charlebois) only work on these changes
and not on the original mVMC package. You can find the original README
and authors of mVMC in the second part of this README file. 

# Citation

You are free to use this code as long as you respect the License terms.
If you use it or learn from it, we ask that you cite the following paper:

[dVMC](https://arxiv.org/abs/1912.09960) or its PRX associated when published.
[mVMC](https://doi.org/10.1016/j.cpc.2018.08.014)

# Usage

Detailed usage is not covered here. Instead many examples can be found
in the "./samples/spectrum/" subdirectory. The easiest way to understand
how to use it is to run these examples. Go there to read the README
and run the few examples. However, the code must be installed first,
see next section.

# Installation:

To install everything, this new version of mVMC must be compiled. 
This can be done using cmake. In the currect (root) directory, do:

$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME
$ make 
$ make install

This should install everything with correct linking in the bin directory 
in your home directory. 

If you have access to intel compilers (icc and ifort), you can impose them 
using the option -DCONFIG=intel in cmake.

It is possible that you do not have the permission to
execute the command "make install" on your cluster.
In that case, make sure you choose an install directory
for which you have permission, in the example above,
we chose $HOME/bin to be the directory of installation.

With this custom path, you need to create the ~/bin/ directory 
and make it execuable by adding this line to the ~/.bashrc 
(if this file does not exist, create it):

export PATH=$PATH:$HOME/bin/

after that, type the command line:

$ source ~/.bashrc 


# Requirements

The libraries required to compile the code are:
- mpi
- lapack or MKL
- blas or MKL
- cmake

It is recommended to consult the original mVMC documentation to learn
more about dependencies and parameters names and idea. Note however that
not all the case in mVMC are covered by the present code, as stated above.


# Details:

Two more README are encouraged to read after this one:

"./samples/spectrum/README" - to learn about the code usage (with few working examples).
"./tool/dvmc/README"        - to learn about the tools that can be used to analyse the data.

Please email me at "maxime.charlebois@usherbrooke.ca" or "maxime.charlebois@uqtr.ca" 
if you have any question on the code. 

This complete the documentation of dVMC.
