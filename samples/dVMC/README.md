<div align="center">

# dVMC calculation sample

</div>

This directory contains a direct application of the dVMC software usage to
compute the ground state of a small 3x4 cluster then the Q-matrix.

## Content

- `params`: Global dVMC solver input parameter file.

- `spectrum_rspace_expected.pdf`: The expected spectrum that any user should
  obtain by running this example.

- `output/`: Working subdirectory.

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

#### Post processing

The expected and relevant output of this sample is the spectrum of the system. To
access it, open the file `spectrum_rspace.pdf`.

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```
