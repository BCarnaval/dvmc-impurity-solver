<div align="center">

# dVMC calculation sample

</div>

This directory contains a direct application of the dVMC software usage to
compute the ground state of a small 3x4 cluster then the Q-matrix.

## Content

- `params`: Global dVMC solver input parameter file.

- `clean.sh`: Cleans this directory except for listed files.

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
dvmc excitations namelist_G.def output/zqp_opt.dat
```

#### Merge binary files

```shell
dvmc process-output output/zvo_nCHAm_nAHCm_0
```

#### Get Q-matrix [^2]

```shell
dvmc qmatrix sqrt
```

[^2]: Add defaults for the three float parameters.
