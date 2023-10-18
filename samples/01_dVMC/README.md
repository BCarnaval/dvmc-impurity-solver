
<div align="center">

# dVMC calculation sample

</div>

This directory contains the minimal code to reproduce the left panel
of Fig. 5 of [arxiv:2209.08092](https://arxiv.org/abs/2209.08092) ([PhysRevB.106.245132](https://doi.org/10.1103/PhysRevB.106.245132)) using a command line interface.

It is a direct application of the dVMC software to compute the 
ground state of a small 2x2 cluster, then the Q-matrix representation
of the Green's function.

## Content

- `params`: Global dVMC solver input parameter file.

- `expected/`: Directory containing expected results from the dVMC calculations.

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

#### Get Q-matrix by computing Green's function

```shell
dvmc green sqrt
```

## Note

You can also use the following command to run all of these in one call

```shell
dvmc run params sqrt
```

and the program will run the 5 steps explained above.

## Post processing

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```
