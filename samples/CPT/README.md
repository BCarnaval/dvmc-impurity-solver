<div align="center">

# CPT calculation sample

</div>

This directory contains a direct application of the dVMC software usage to
compute the ground state of a small 3x4 cluster then using the Q-matrix
representation in cluster perturbation theory (CPT) to output the Fermi
surface of the system by calling PyQCM library.

## Content

- `model_3x4.py`: Defines the lattice model using PyQCM library instances.

- `fermi_surface.py`: Reads the Q-matrix representation of the ground state
  to generate the Fermi surface using PyQCM library.

- `params`: Global dVMC solver input parameter file.

- `QCM_param.def`[^1]: PyQCM specificly formatted parameter file.

- `fermi_surface_expected.pdf`: The expected Fermi surface that any user should
  obtain by running this example.

- `output/`: Working subdirectory.

[^1]: Might be yeeted soon, it's content is contained inside params file.

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
