
<div align="center">

# CPT-dVMC calculation sample

</div>

This directory contains the minimal code to reproduce the results of
Fig. 3a) and the third pannel of Fig. 5 of [arxiv:2209.08092](https://arxiv.org/abs/2209.08092) 
([PhysRevB.106.245132](https://doi.org/10.1103/PhysRevB.106.245132)) using PyQCM.

It computes the ground state, density of states and Fermi surface of 
a small 3x4 cluster then using the Q-matrix representation in cluster 
perturbation theory (CPT) to output the cluster spectral function 
of the system by calling PyQCM library.

> ### Note
>
> It takes approximately 40 minutes to run [Béluga](https://docs.alliancecan.ca/wiki/B%C3%A9luga/en)
> supercomputer using 64 MPI processes with 1 CPU per task and 1G of memory. It roughly reproduces the
>  third graph of figure 3 in the paper: **[PhysRevB.106.245132](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.245132)**

## Content

- `model_3x4.py`: Defines the lattice model using PyQCM library instances.

- `cluster_spectral.py`: Uses PyQCM dVMC interface to find the groundstate then
  CPT method to output the cluster spectral function of the system.

- `params`: Global dVMC solver input parameter file.

- `expected/`: Directory containing expected results from the CPT-dVMC calculations.

## Usage

Since this example uses dVMC through PyQCM, the only thing to do is to execute
the `cluster_spectral.py` script using Python 3 (with proper PyQCM installation)

#### Run script to generate cluster spectral function from CPT-dVMC

```shell
python3 cluster_spectral.py
```

## Post processing

To clean the directory from generated files that would be overwritten by the program,
use the `dvmc` command line interface

```shell
dvmc clean --help
```
