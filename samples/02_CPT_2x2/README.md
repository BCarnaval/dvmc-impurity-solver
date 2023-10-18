
<div align="center">

# CPT-dVMC calculation sample

</div>

This directory contains the minimal code to reproduce the left panel
of Fig. 5 of [arxiv:2209.08092](https://arxiv.org/abs/2209.08092) ([PhysRevB.106.245132](https://doi.org/10.1103/PhysRevB.106.245132)) using PyQCM.

It is a direct application of the dVMC software to compute the ground 
state of a small 2x2 cluster, then using the Q-matrix representation in 
cluster perturbation theory (CPT) to output the cluster spectral function, 
fermi surface and dos of the system by calling PyQCM.

This model is the same as the example contained in sample `01_dVMC`, but with
more options due to the PyQCM machinery (fermi_surface and spectrum).

## Content

- `model_2x2.py`: Defines the lattice model using PyQCM library instances.

- `cluster_spectral.py`: Uses PyQCM dVMC interface to find the groundstate then
  CPT method to output the cluster spectral function.

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
