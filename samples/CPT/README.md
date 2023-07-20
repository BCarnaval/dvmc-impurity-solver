<div align="center">

# CPT calculation sample

</div>

This directory contains a direct application of the dVMC software to
compute the ground state of a small 3x4 cluster then using the Q-matrix
representation in cluster perturbation theory (CPT) to output the Fermi
surface of the system by calling PyQCM library.

## Content

- `model_3x4.py`: Defines the lattice model using PyQCM library instances.

- `fermi_surface.py`: Reads the Q-matrix representation of the ground state
  to generate the Fermi surface using PyQCM library.

- `params`: Global dVMC solver input parameter file.

- `QCM_param.def`: PyQCM specificly formatted parameter file.

- `expected/`: Directory containing expected results from the CPT calculations.

## Usage

Here a list of commands to use considering that the CMake installation has been
successful to run the sample:

#### Generate input files

```shell
dvmc generate params 3
```

#### Run dVMC to get ground state

```shell
dvmc groundstate
```

#### Run dVMC to get excitations

```shell
dvmc excitations
```

#### Merge binary files

```shell
dvmc process-output output/zvo_nCHAm_nAHCm_0
```

#### Get Q-matrix

```shell
dvmc qmatrix sqrt
```

#### Run script to generate A(k,w) from CPT

```shell
python3 fermi_surface.py output/qmatrix.def
```

## Post processing

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```
